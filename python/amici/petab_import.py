"""
Import a model in the PEtab (https://github.com/ICB-DCM/PEtab/) format into
AMICI.
"""

import argparse
import logging
import math
import os
import time
from typing import List, Dict, Union, Optional

import amici
import libsbml
import pandas as pd
import petab
import sympy as sp
from amici.logging import get_logger
from petab import (OBSERVABLE_TRANSFORMATION, NOISE_DISTRIBUTION,
                   CONDITION_NAME, OBSERVABLE_ID)


logger = get_logger(__name__, logging.WARNING)


def get_fixed_parameters(
        sbml_model: 'libsbml.Model',
        condition_df: Optional[pd.DataFrame] = None,
        const_species_to_parameters: bool = False) -> List[str]:
    """Determine, set and return fixed model parameters

    Parameters specified in `condition_df` are turned into constants.
    Only global SBML parameters are considered. Local parameters are ignored.

    Arguments:
        condition_df:
            PEtab condition table. If provided, the respective parameters
            will be turned into AMICI constant parameters.
        sbml_model:
            libsbml.Model instance
        const_species_to_parameters:
            If `True`, species which are marked constant within the SBML model
            will be turned into constant parameters *within* the given
            `sbml_model`.

    Returns:
        List of IDs of parameters which are to be considered constant.
    """

    # Column names are model parameter IDs, compartment IDs or species IDs.
    # Thereof, all parameters except for any overridden ones should be made
    # constant.
    # (Could potentially still be made constant, but leaving them might
    # increase model reusability)

    # handle parameters in condition table
    if condition_df is not None:
        fixed_parameters = list(condition_df.columns)
        # get rid of conditionName column
        try:
            fixed_parameters.remove(CONDITION_NAME)
        except ValueError:
            pass

        # remove overridden parameters (`object`-type columns)
        fixed_parameters = [p for p in fixed_parameters
                            if condition_df[p].dtype != 'O'
                            and sbml_model.getParameter(p) is not None]
        # must be unique
        if len(fixed_parameters) != len(set(fixed_parameters)):
            raise AssertionError(
                'len(fixed_parameters) != len(set(fixed_parameters))')
    else:
        fixed_parameters = []

    # Others are optional
    if const_species_to_parameters:
        # Turn species which are marked constant in the SBML model into
        # parameters
        constant_species = constant_species_to_parameters(sbml_model)

        logger.debug("Constant species converted to parameters: "
                     + str(len(constant_species)))
        logger.info("Non-constant species "
                    + str(len(sbml_model.getListOfSpecies())))

        # ... and append them to the list of fixed_parameters
        for species in constant_species:
            if species not in fixed_parameters:
                fixed_parameters.append(species)

    # Ensure mentioned parameters exist in the model. Remove additional ones
    # from list
    for fixed_parameter in fixed_parameters[:]:
        # check global parameters
        if not sbml_model.getParameter(fixed_parameter) \
                and not sbml_model.getSpecies(fixed_parameter):
            logger.warning(f"Parameter or species '{fixed_parameter}'"
                           " provided in condition table but not present in"
                           " model. Ignoring.")
            fixed_parameters.remove(fixed_parameter)

    if condition_df is None:
        return fixed_parameters

    # initial concentrations for species or initial compartment sizes in
    # condition table will need to be turned into fixed parameters

    # if there is no initial assignment for that species, we'd need
    # to create one. to avoid any naming collision right away, we don't allow
    # that for now

    # we can't handle them yet
    compartments = [col for col in condition_df
               if sbml_model.getCompartment(col) is not None]
    if compartments:
        raise NotImplementedError("Can't handle initial compartment sizes "
                                  "at the moment. Consider creating an "
                                  f"initial assignment for {compartments}")

    species = [col for col in condition_df
               if sbml_model.getSpecies(col) is not None]
    if species:
        raise NotImplementedError("Can't handle species in condition table."
                                  "Consider creating an initial assignment for"
                                  f" {species}")

    return fixed_parameters


def species_to_parameters(species_ids: List[str],
                          sbml_model: 'libsbml.Model') -> List[str]:
    """Turn a SBML species into parameters and replace species references
    inside the model instance.

    Arguments:
        species_ids: List of SBML species ID to convert to parameters with the
            same ID as the replaced species.
        sbml_model: SBML model to modify

    Returns:
        List of IDs of species which have been converted to parameters
    """
    transformables = []

    for species_id in species_ids:
        species = sbml_model.getSpecies(species_id)

        if species.getHasOnlySubstanceUnits():
            logger.warning(
                f"Ignoring {species.getId()} which has only substance units."
                " Conversion not yet implemented.")
            continue

        if math.isnan(species.getInitialConcentration()):
            logger.warning(
                f"Ignoring {species.getId()} which has no initial "
                "concentration. Amount conversion not yet implemented.")
            continue

        transformables.append(species_id)

    # Must not remove species while iterating over getListOfSpecies()
    for species_id in transformables:
        species = sbml_model.removeSpecies(species_id)
        par = sbml_model.createParameter()
        par.setId(species.getId())
        par.setName(species.getName())
        par.setConstant(True)
        par.setValue(species.getInitialConcentration())
        par.setUnits(species.getUnits())

    # Remove from reactants and products
    for reaction in sbml_model.getListOfReactions():
        for species_id in transformables:
            # loop, since removeX only removes one instance
            while reaction.removeReactant(species_id):
                pass
            while reaction.removeProduct(species_id):
                pass
            while reaction.removeModifier(species_id):
                pass

    return transformables


def constant_species_to_parameters(sbml_model: 'libsbml.Model') -> List[str]:
    """Convert constant species in the SBML model to constant parameters.

    This can be used e.g. for setting up models with condition-specific
    constant species for PEtab, since there it is not possible to specify
    constant species in the condition table.

    Arguments:
        sbml_model: SBML Model

    Returns:
        List of IDs of SBML species that have been turned into constants

    Raises:

    """
    transformables = []
    for species in sbml_model.getListOfSpecies():
        if not species.getConstant() and not species.getBoundaryCondition():
            continue

        transformables.append(species.getId())

    return species_to_parameters(transformables, sbml_model)


def import_model(sbml_model: Union[str, 'libsbml.Model'],
                 condition_table: Optional[Union[str, pd.DataFrame]] = None,
                 measurement_table: Optional[Union[str, pd.DataFrame]] = None,
                 model_name: Optional[str] = None,
                 model_output_dir: Optional[str] = None,
                 verbose: bool = True,
                 allow_reinit_fixpar_initcond: bool = True,
                 **kwargs) -> None:
    """Create AMICI model from PEtab problem

    Arguments:
        sbml_model:
            PEtab SBML model or SBML file name.
        condition_table:
            PEtab condition table. If provided, parameters from there will be
            turned into AMICI constant parameters (i.e. parameters w.r.t. which
            no sensitivities will be computed).
        measurement_table:
            PEtab measurement table.
        model_name:
            Name of the generated model. If model file name was provided,
            this defaults to the file name without extension, otherwise
            the SBML model ID will be used.
        model_output_dir:
            Directory to write the model code to. Will be created if doesn't
            exist. Defaults to current directory.
        verbose:
            Print/log extra information.
        allow_reinit_fixpar_initcond:
            See amici.ode_export.ODEExporter. Must be enabled if initial
            states are to be reset after preequilibration.
        **kwargs:
            Additional keyword arguments to be passed to
            ``amici.sbml_importer.sbml2amici``.
    """
    if verbose:
        logger.setLevel(verbose)

    logger.info(f"Importing model ...")

    # Create a copy, because it will be modified by SbmlImporter
    sbml_doc = sbml_model.getSBMLDocument().clone()
    sbml_model = sbml_doc.getModel()

    sbml_importer = amici.SbmlImporter(sbml_model)

    # Model name from SBML ID or filename
    if model_name is None:
        if isinstance(sbml_model, libsbml.Model):
            model_name = sbml_model.getId()
        else:
            model_name = os.path.splitext(os.path.split(sbml_model)[-1])[0]

    if model_output_dir is None:
        model_output_dir = os.path.join(os.getcwd(), model_name)

    sbml_model = sbml_importer.sbml

    logger.info(f"Model name is '{model_name}'. "
                f"Writing model code to '{model_output_dir}'.")
    show_model_info(sbml_model)

    # Read PEtab observables and sigmas
    observables = petab.get_observables(sbml_model, remove=True)
    sigmas = petab.get_sigmas(sbml_model, remove=True)

    # Read PEtab error model
    if measurement_table is not None:
        if isinstance(measurement_table, str):
            measurement_df = petab.get_measurement_df(measurement_table)
        else:
            measurement_df = measurement_table

        noise_distrs = petab_noise_distributions_to_amici(
            petab.get_noise_distributions(measurement_df))
    else:
        # use default
        noise_distrs = {}

    # Replace observableIds occurring in error model definition
    for observable_id, formula in sigmas.items():
        repl = sp.sympify(formula).subs(
            observable_id, observables[observable_id]['formula'])
        sigmas[observable_id] = str(repl)

    logger.info(f'Observables: {len(observables)}')
    logger.info(f'Sigmas: {len(sigmas)}')

    if not len(sigmas) == len(observables):
        raise AssertionError(
            f'Number of provided observables ({len(observables)}) and sigmas '
            f'({len(sigmas)}) do not match.')

    if condition_table is not None:
        # get the condition dataframe before parsing fixed parameters
        if isinstance(condition_table, str):
            condition_df = petab.get_condition_df(condition_table)
        else:
            condition_df = condition_table
        logger.debug(f'Condition table: {condition_df.shape}')
    else:
        condition_df = None

    fixed_parameters = get_fixed_parameters(sbml_model=sbml_model,
                                            condition_df=condition_df)

    logger.debug(f"Fixed parameters are {fixed_parameters}")
    logger.info(f"Overall fixed parameters: {len(fixed_parameters)}")
    logger.info("Non-constant global parameters: "
                + str(len(sbml_model.getListOfParameters())
                      - len(fixed_parameters)))

    # Create Python module from SBML model
    start = time.time()
    sbml_importer.sbml2amici(
        modelName=model_name,
        output_dir=model_output_dir,
        observables=observables,
        constantParameters=fixed_parameters,
        sigmas=sigmas,
        allow_reinit_fixpar_initcond=allow_reinit_fixpar_initcond,
        noise_distributions=noise_distrs,
        verbose=verbose,
        **kwargs)
    end = time. time()

    logger.info(f"Model imported successfully in {round(end - start, 2)}s")


def petab_noise_distributions_to_amici(noise_distributions: Dict) -> Dict:
    """
    Map from the petab to the amici format of noise distribution
    identifiers.

    Arguments:
        noise_distributions: as obtained from `petab.get_noise_distributions`

    Returns:
        Dictionary of observable_id => AMICI noise-distributions
    """
    amici_distrs = {}
    for id_, val in noise_distributions.items():
        amici_val = ''

        if val[OBSERVABLE_TRANSFORMATION]:
            amici_val += val[OBSERVABLE_TRANSFORMATION] + '-'

        if val[NOISE_DISTRIBUTION]:
            amici_val += val[NOISE_DISTRIBUTION]

        amici_distrs[id_] = amici_val

    return amici_distrs


def petab_scale_to_amici_scale(scale_str: str) -> int:
    """Convert PEtab parameter scaling string to AMICI scaling integer"""

    if scale_str == petab.LIN:
        return amici.ParameterScaling_none
    if scale_str == petab.LOG:
        return amici.ParameterScaling_ln
    if scale_str == petab.LOG10:
        return amici.ParameterScaling_log10

    raise ValueError(f"Invalid parameter scale {scale_str}")


def show_model_info(sbml_model: 'libsbml.Model'):
    """Log some model quantities"""

    logger.info(f'Species: {len(sbml_model.getListOfSpecies())}')
    logger.info('Global parameters: '
                + str(len(sbml_model.getListOfParameters())))
    logger.info(f'Reactions: {len(sbml_model.getListOfReactions())}')


def parse_cli_args():
    """Parse command line arguments

    Returns:
        Parsed CLI arguments from ``argparse``.
    """

    parser = argparse.ArgumentParser(
        description='Import PEtab-format model into AMICI.')

    # General options:
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='More verbose output')
    parser.add_argument('-o', '--output-dir', dest='model_output_dir',
                        help='Name of the model directory to create')
    parser.add_argument('--no-compile', action='store_false',
                        dest='compile',
                        help='Only generate model code, do not compile')

    # Call with set of files
    parser.add_argument('-s', '--sbml', dest='sbml_file_name',
                        help='SBML model filename')
    parser.add_argument('-m', '--measurements', dest='measurement_file_name',
                        help='Measurement table')
    parser.add_argument('-c', '--conditions', dest='condition_file_name',
                        help='Conditions table')
    parser.add_argument('-p', '--parameters', dest='parameter_file_name',
                        help='Parameter table')

    parser.add_argument('-y', '--yaml', dest='yaml_file_name',
                       help='PEtab YAML problem filename')

    parser.add_argument('-n', '--model-name', dest='model_name',
                        help='Name of the python module generated for the '
                             'model')

    args = parser.parse_args()

    if not args.yaml_file_name \
            and not all((args.sbml_file_name, args.condition_file_name,
                         args.measurement_file_name)):
        parser.error('When not specifying a model name or YAML file, then '
                     'SBML, condition and measurement file must be specified')

    return args


def main():
    """
    Command line interface to import a model in the PEtab
    (https://github.com/ICB-DCM/PEtab/) format into AMICI.
    """
    args = parse_cli_args()

    if args.yaml_file_name:
        pp = petab.Problem.from_yaml(args.yaml_file_name)
    else:
        pp = petab.Problem.from_files(
            sbml_file=args.sbml_file_name,
            condition_file=args.condition_file_name,
            measurement_file=args.measurement_file_name,
            parameter_file=args.parameter_file_name)

    # First check for valid PEtab
    petab.lint_problem(pp)

    import_model(model_name=args.model_name,
                 sbml_model=pp.sbml_model,
                 condition_table=pp.condition_df,
                 measurement_table=pp.measurement_df,
                 model_output_dir=args.model_output_dir,
                 compile=args.compile,
                 verbose=args.verbose)


if __name__ == '__main__':
    main()
