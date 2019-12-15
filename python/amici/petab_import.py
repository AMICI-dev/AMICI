"""
Import a model in the PEtab (https://github.com/ICB-DCM/PEtab/) format into
AMICI
"""

import amici
import os
import time
import math
import logging
from typing import List, Dict, Union, Optional
import pandas as pd
import argparse

import libsbml
import petab
from colorama import Fore
from colorama import init as init_colorama


logger = logging.getLogger(__name__)


def get_fixed_parameters(
        sbml_model: 'libsbml.Model',
        condition_df: Optional[pd.DataFrame] = None,
        const_species_to_parameters: bool = False) -> List[str]:
    """Determine, set and return fixed model parameters

    Parameters specified in `condition_file_name` are turned into constants.
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
        List of IDs of parameter which are to be considered constant
    """

    # column names are model parameter names that should be made constant
    # except for any overridden parameters
    # (Could potentially still be made constant, but leaving them might
    # increase model reusability)
    if condition_df is not None:
        fixed_parameters = list(condition_df.columns)
        try:
            fixed_parameters.remove('conditionName')
        except ValueError:
            pass
        # remove overridden parameters
        fixed_parameters = [p for p in fixed_parameters
                            if condition_df[p].dtype != 'O']
        # must be unique
        assert(len(fixed_parameters) == len(set(fixed_parameters)))
    else:
        fixed_parameters = []

    # States occurring as column names of the condition table need to be
    #  converted to parameters
    # TODO https://github.com/ICB-DCM/PEtab/issues/181
    species_to_convert = [x for x in fixed_parameters
                          if sbml_model.getSpecies(x)]
    species_to_parameters(species_to_convert, sbml_model)

    # Others are optional
    if const_species_to_parameters:
        # Turn species which are marked constant in the SBML model into
        # parameters
        constant_species = constant_species_to_parameters(sbml_model)

        logger.log(logging.INFO, "Constant species converted to parameters "
                   + str(len(constant_species)))
        logger.log(logging.INFO, "Non-constant species "
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
            logger.log(logging.WARN,
                       f"{Fore.YELLOW}Parameter or species '{fixed_parameter}'"
                       " provided in condition table but not present in"
                       " model.")
            fixed_parameters.remove(fixed_parameter)

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
            return

        if math.isnan(species.getInitialConcentration()):
            logger.warning(
                f"Ignoring {species.getId()} which has no initial "
                "concentration. Amount conversion not yet implemented.")
            return

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
                 model_output_dir: str = None,
                 verbose: bool = True,
                 allow_reinit_fixpar_initcond: bool = False,
                 **kwargs) -> None:
    """Create AMICI model from PEtab problem

    Arguments:
        sbml_model:
            PEtab SBML model. If `libsbml.Model` is passed, this will be
            modified inplace.
        condition_table:
            PEtab condition table. If provided, parameters from there will be
            turned into AMICI constant parameters (i.e. parameter w.r.t. which
            no sensitivities will be computed).
        measurement_table:
            PEtab measurement table
        model_name:
            Name of the generated model. If model file name was provided,
            this defaults to the file name without extension, otherwise
            the SBML model ID will be used.
        model_output_dir:
            Directory to write the model code to. Will be created if doesn't
            exist. Defaults to current directory.
        verbose:
            Print/log extra information
        allow_reinit_fixpar_initcond:
            see amici.ode_export.ODEExporter
        **kwargs:
            Additional keyword arguments to be passed to
            ``amici.sbml_importer.sbml2amici``.
    """

    # Color output
    init_colorama(autoreset=True)

    if verbose:
        logger.log(logging.INFO, f"{Fore.GREEN}Importing model ...")

    sbml_importer = amici.SbmlImporter(sbml_model)

    if model_name is None:
        if isinstance(sbml_model, libsbml.Model):
            model_name = sbml_model.getId()
        else:
            model_name = os.path.splitext(os.path.split(sbml_model)[-1])[0]

    if model_output_dir is None:
        model_output_dir = os.path.join(os.getcwd(), model_name)

    sbml_model = sbml_importer.sbml

    if verbose:
        logger.log(logging.INFO,
                   f"{Fore.GREEN}Model name is '{model_name}' "
                   f"Writing model code to '{model_output_dir}'")
        show_model_info(sbml_model)

    # Read PEtab observables and sigmas
    observables = petab.get_observables(sbml_importer.sbml, remove=True)
    sigmas = petab.get_sigmas(sbml_importer.sbml, remove=True)

    # Read PEtab error model
    if measurement_table is not None:
        if isinstance(measurement_table, str):
            measurement_df = petab.get_measurement_df(measurement_table)
        else:
            measurement_df = measurement_table

        noise_distrs = petab_noise_distributions_to_amici(
            petab.get_noise_distributions(measurement_df))
    else:
        noise_distrs = {}  # use default

    # Replace observablesIds occurring in error model definition
    import sympy as sp
    for observable_id, formula in sigmas.items():
        repl = sp.sympify(formula).subs(
            observable_id, observables[observable_id]['formula'])
        sigmas[observable_id] = str(repl)

    if verbose:
        logger.log(logging.INFO, f'Observables {len(observables)}')
        logger.log(logging.INFO, f'Sigmas {len(sigmas)}')

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
        logger.log(logging.INFO, f'Condition table: {condition_df.shape}')
    else:
        condition_df = None

    fixed_parameters = get_fixed_parameters(sbml_model=sbml_model,
                                            condition_df=condition_df)

    if verbose:
        logger.log(logging.INFO,
                   f"Overall fixed parameters {len(fixed_parameters)}")
        logger.log(logging.INFO, "Non-constant global parameters "
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
        **kwargs)
    end = time. time()

    if verbose:
        logger.log(logging.INFO, f"{Fore.GREEN}Model imported successfully in "
                                 f"{round(end - start, 2)}s")


def petab_noise_distributions_to_amici(noise_distributions: Dict) -> Dict:
    """
    Map from the petab to the amici format of noise distribution
    identifiers.

    Arguments:
        noise_distributions: as obtained from `petab.get_noise_distributions`

    Returns:
        Dictionary of obserable_id => AMICI noise-distributions
    """
    amici_distrs = {}
    for id_, val in noise_distributions.items():
        amici_val = ''

        if val['observableTransformation']:
            amici_val += val['observableTransformation'] + '-'

        if val['noiseDistribution']:
            amici_val += val['noiseDistribution']

        amici_distrs[id_] = amici_val

    return amici_distrs


def show_model_info(sbml_model: 'libsbml.Model'):
    """Log some model quantities"""

    logger.log(logging.INFO, f'Species: {len(sbml_model.getListOfSpecies())}')
    logger.log(logging.INFO, 'Global parameters: '
               + str(len(sbml_model.getListOfParameters())))
    logger.log(logging.INFO,
               f'Reactions: {len(sbml_model.getListOfReactions())}')


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

    # or with model name, following default naming
    parser.add_argument('-n', '--model-name', dest='model_name',
                        help='Model name where all files are in the working '
                        'directory and follow PEtab naming convention. '
                        'Specifying -[smcp] will override defaults')

    args = parser.parse_args()

    if args.model_name:
        if not args.sbml_file_name:
            args.sbml_file_name = petab.get_default_sbml_file_name(
                args.model_name)
        if not args.measurement_file_name:
            args.measurement_file_name = \
                petab.get_default_measurement_file_name(args.model_name)
        if not args.condition_file_name:
            args.condition_file_name = petab.get_default_condition_file_name(
                args.model_name)
        if not args.parameter_file_name:
            args.parameter_file_name = petab.get_default_parameter_file_name(
                args.model_name)

    if not args.model_name and \
            (not args.sbml_file_name
             or not args.condition_file_name
             or not args.measurement_file_name):
        parser.error('When not specifying a model name, sbml, '
                     'condition and measurement file must be specified')

    return args


def main():
    """
    Command line interface to import a model in the PEtab
    (https://github.com/ICB-DCM/PEtab/) format into AMICI.
    """
    args = parse_cli_args()

    # First check for valid PEtab
    pp = petab.Problem.from_files(
        sbml_file=args.sbml_file_name,
        condition_file=args.condition_file_name,
        measurement_file=args.measurement_file_name,
        parameter_file=args.parameter_file_name)
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
