"""
PEtab Import
------------
Import a model in the :mod:`petab` (https://github.com/ICB-DCM/PEtab/) format
into AMICI.
"""

import argparse
import importlib
import logging
import math
import os
import shutil
import sys
import tempfile
from typing import List, Dict, Union, Optional, Tuple

import amici
import libsbml
import numpy as np
import pandas as pd
import petab
import sympy as sp
from amici.logging import get_logger, log_execution_time, set_log_level
from petab.C import *

logger = get_logger(__name__, logging.WARNING)


def get_fixed_parameters(
        sbml_model: 'libsbml.Model',
        condition_df: Optional[pd.DataFrame] = None,
        const_species_to_parameters: bool = False) -> List[str]:
    """
    Determine, set and return fixed model parameters.

    Parameters specified in `condition_df` are turned into constants.
    Only global SBML parameters are considered. Local parameters are ignored.

    :param condition_df:
        PEtab condition table. If provided, the respective parameters
        will be turned into AMICI constant parameters.

    :param sbml_model:
        libsbml.Model instance

    :param const_species_to_parameters:
        If `True`, species which are marked constant within the SBML model
        will be turned into constant parameters *within* the given
        `sbml_model`.

    :return:
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

        logger.debug(f'Condition table: {condition_df.shape}')

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
               if not np.issubdtype(condition_df[col].dtype, np.number)
               and sbml_model.getSpecies(col) is not None]
    if species:
        raise NotImplementedError(
            "Can't handle parameterized initial concentrations in condition "
            f"table. Consider creating an initial assignment for {species}")

    return fixed_parameters


def species_to_parameters(species_ids: List[str],
                          sbml_model: 'libsbml.Model') -> List[str]:
    """
    Turn a SBML species into parameters and replace species references
    inside the model instance.

    :param species_ids:
        List of SBML species ID to convert to parameters with the same ID as
        the replaced species.

    :param sbml_model:
        SBML model to modify

    :return:
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
    """
    Convert constant species in the SBML model to constant parameters.

    This can be used e.g. for setting up models with condition-specific
    constant species for PEtab, since there it is not possible to specify
    constant species in the condition table.

    :param sbml_model:
        SBML Model

    :return:
        List of IDs of SBML species that have been turned into constants
    """
    transformables = []
    for species in sbml_model.getListOfSpecies():
        if not species.getConstant() and not species.getBoundaryCondition():
            continue

        transformables.append(species.getId())

    return species_to_parameters(transformables, sbml_model)


def import_petab_problem(
        petab_problem: petab.Problem,
        model_output_dir: str = None,
        model_name: str = None,
        force_compile: bool = False,
        **kwargs) -> 'amici.Model':
    """
    Import model from petab problem.

    :param petab_problem:
        A petab problem containing all relevant information on the model.

    :param model_output_dir:
        Directory to write the model code to. Will be created if doesn't
        exist. Defaults to current directory.

    :param model_name:
        Name of the generated model. If model file name was provided,
        this defaults to the file name without extension, otherwise
        the SBML model ID will be used.

    :param force_compile:
        Whether to compile the model even if the target folder is not empty,
        or the model exists already.

    :param kwargs:
        Additional keyword arguments to be passed to
        :meth:`amici.sbml_import.SbmlImporter.sbml2amici`.

    :return:
        The imported model.
    """
    # generate folder and model name if necessary
    if model_output_dir is None:
        model_output_dir = \
            _create_model_output_dir_name(petab_problem.sbml_model)
    else:
        model_output_dir = os.path.abspath(model_output_dir)

    if model_name is None:
        model_name = _create_model_name(model_output_dir)

    # create folder
    if not os.path.exists(model_output_dir):
        os.makedirs(model_output_dir)

    # add to path
    if model_output_dir not in sys.path:
        sys.path.insert(0, model_output_dir)

    # check if compilation necessary
    if not _can_import_model(model_name) or force_compile:
        # check if folder exists
        if os.listdir(model_output_dir) and not force_compile:
            raise ValueError(
                f"Cannot compile to {model_output_dir}: not empty. "
                "Please assign a different target or set `force_compile`.")

        # remove folder if exists
        if os.path.exists(model_output_dir):
            shutil.rmtree(model_output_dir)

        logger.info(f"Compiling model {model_name} to {model_output_dir}.")

        # compile the model
        import_model(sbml_model=petab_problem.sbml_model,
                     condition_table=petab_problem.condition_df,
                     observable_table=petab_problem.observable_df,
                     model_name=model_name,
                     model_output_dir=model_output_dir,
                     **kwargs)
        # ensure we will find the newly created module
        importlib.invalidate_caches()

    # load module
    if model_name in sys.modules:
        # reload, because may just have been created
        importlib.reload(sys.modules[model_name])
        model_module = sys.modules[model_name]
    else:
        model_module = importlib.import_module(model_name)

    # import model
    model = model_module.getModel()

    logger.info(f"Successfully loaded model {model_name} "
                f"from {model_output_dir}.")

    return model


def _create_model_output_dir_name(sbml_model: 'libsbml.Model') -> str:
    """
    Find a folder for storing the compiled amici model.
    If possible, use the sbml model id, otherwise create a random folder.
    The folder will be located in the `amici_models` subfolder of the current
    folder.
    """
    BASE_DIR = os.path.abspath("amici_models")

    # create base directory
    if not os.path.exists(BASE_DIR):
        os.makedirs(BASE_DIR)

    # try sbml model id
    sbml_model_id = sbml_model.getId()
    if sbml_model_id:
        model_output_dir = os.path.join(BASE_DIR, sbml_model_id)
    else:
        # create random folder name
        model_output_dir = tempfile.mkdtemp(dir=BASE_DIR)

    return model_output_dir


def _create_model_name(folder: str) -> str:
    """
    Create a name for the model.
    Just re-use the last part of the folder.
    """
    return os.path.split(os.path.normpath(folder))[-1]


def _can_import_model(model_name: str) -> bool:
    """
    Check whether a module of that name can already be imported.
    """
    # try to import (in particular checks version)
    try:
        importlib.import_module(model_name)
    except ModuleNotFoundError:
        return False

    # no need to (re-)compile
    return True


@log_execution_time('Importing PEtab model', logger)
def import_model(sbml_model: Union[str, 'libsbml.Model'],
                 condition_table: Optional[Union[str, pd.DataFrame]] = None,
                 observable_table: Optional[Union[str, pd.DataFrame]] = None,
                 model_name: Optional[str] = None,
                 model_output_dir: Optional[str] = None,
                 verbose: Optional[Union[bool,int]] = True,
                 allow_reinit_fixpar_initcond: bool = True,
                 **kwargs) -> None:
    """
    Create AMICI model from PEtab problem

    :param sbml_model:
        PEtab SBML model or SBML file name.

    :param condition_table:
        PEtab condition table. If provided, parameters from there will be
        turned into AMICI constant parameters (i.e. parameters w.r.t. which
        no sensitivities will be computed).

    :param observable_table:
        PEtab observable table.

    :param model_name:
        Name of the generated model. If model file name was provided,
        this defaults to the file name without extension, otherwise
        the SBML model ID will be used.

    :param model_output_dir:
        Directory to write the model code to. Will be created if doesn't
        exist. Defaults to current directory.

    :param verbose:
        Print/log extra information.

    :param allow_reinit_fixpar_initcond:
        See :class:`amici.ode_export.ODEExporter`. Must be enabled if initial
        states are to be reset after preequilibration.

    :param kwargs:
        Additional keyword arguments to be passed to
        :meth:`amici.sbml_import.SbmlImporter.sbml2amici`.
    """

    set_log_level(logger, verbose)

    logger.info(f"Importing model ...")

    # Get PEtab tables
    observable_df = petab.get_observable_df(observable_table)
    # to determine fixed parameters
    condition_df = petab.get_condition_df(condition_table)

    if observable_df is None:
        raise NotImplementedError("PEtab import without observables table "
                                  "is currently not supported.")

    # Model name from SBML ID or filename
    if model_name is None:
        if isinstance(sbml_model, libsbml.Model):
            model_name = sbml_model.getId()
        else:
            model_name = os.path.splitext(os.path.split(sbml_model)[-1])[0]

    if model_output_dir is None:
        model_output_dir = os.path.join(os.getcwd(), model_name)

    logger.info(f"Model name is '{model_name}'. "
                f"Writing model code to '{model_output_dir}'.")

    # Load model
    if isinstance(sbml_model, str):
        # from file
        sbml_reader = libsbml.SBMLReader()
        sbml_doc = sbml_reader.readSBMLFromFile(sbml_model)
        sbml_model = sbml_doc.getModel()
    else:
        # Create a copy, because it will be modified by SbmlImporter
        sbml_doc = sbml_model.getSBMLDocument().clone()
        sbml_model = sbml_doc.getModel()

    show_model_info(sbml_model)

    sbml_importer = amici.SbmlImporter(sbml_model)
    sbml_model = sbml_importer.sbml

    if observable_df is not None:
        observables, noise_distrs, sigmas = \
            get_observation_model(observable_df)

    logger.info(f'Observables: {len(observables)}')
    logger.info(f'Sigmas: {len(sigmas)}')

    if not len(sigmas) == len(observables):
        raise AssertionError(
            f'Number of provided observables ({len(observables)}) and sigmas '
            f'({len(sigmas)}) do not match.')

    # TODO: adding extra output parameters is currently not supported,
    #  so we add any output parameters to the SBML model.
    #  this should be changed to something more elegant
    # <BeginWorkAround>
    formulas = {val['formula'] for val in observables.values()}
    formulas |= set(sigmas.values())
    output_parameters = set()
    for formula in formulas:
        for free_sym in sp.sympify(formula).free_symbols:
            sym = str(free_sym)
            if sbml_model.getElementBySId(sym) is None:
                output_parameters.add(sym)
    logger.debug(f"Adding output parameters to model: {output_parameters}")
    for par in output_parameters:
        petab.add_global_parameter(sbml_model, par)
    # <EndWorkAround>

    fixed_parameters = get_fixed_parameters(sbml_model=sbml_model,
                                            condition_df=condition_df)

    logger.debug(f"Fixed parameters are {fixed_parameters}")
    logger.info(f"Overall fixed parameters: {len(fixed_parameters)}")
    logger.info("Variable parameters: "
                + str(len(sbml_model.getListOfParameters())
                      - len(fixed_parameters)))

    # Create Python module from SBML model
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


def get_observation_model(observable_df: pd.DataFrame
                          ) -> Tuple[Dict[str, Dict[str, str]],
                                     Dict[str, str],
                                     Dict[str, Union[str, float]]]:
    """
    Get observables, sigmas, and noise distributions from PEtab observation
    table in a format suitable for
    :meth:`amici.sbml_import.SbmlImporter.sbml2amici`.

    :param observable_df:
        PEtab observables table

    :return:
        Tuple of dicts with observables, noise distributions, and sigmas.
    """

    if observable_df is None:
        return {}, {}, {}

    observables = {}
    sigmas = {}

    for _, observable in observable_df.iterrows():
        oid = observable.name
        name = observable.get(OBSERVABLE_NAME, "")
        formula_obs = observable[OBSERVABLE_FORMULA]
        formula_noise = observable[NOISE_FORMULA]
        observables[oid] = {'name': name, 'formula': formula_obs}
        sigmas[oid] = formula_noise

    # Replace observableIds occurring in error model definition
    for observable_id, formula in sigmas.items():
        repl = sp.sympify(formula).subs(
            observable_id, observables[observable_id]['formula'])
        sigmas[observable_id] = str(repl)

    noise_distrs = petab_noise_distributions_to_amici(observable_df)

    return observables, noise_distrs, sigmas


def petab_noise_distributions_to_amici(observable_df: pd.DataFrame) -> Dict:
    """
    Map from the petab to the amici format of noise distribution
    identifiers.

    :param observable_df:
        PEtab observable table

    :return:
        Dictionary of observable_id => AMICI noise-distributions
    """
    amici_distrs = {}
    for _, observable in observable_df.iterrows():
        amici_val = ''

        if OBSERVABLE_TRANSFORMATION in observable \
                and isinstance(observable[OBSERVABLE_TRANSFORMATION], str) \
                and observable[OBSERVABLE_TRANSFORMATION]:
            amici_val += observable[OBSERVABLE_TRANSFORMATION] + '-'

        if NOISE_DISTRIBUTION in observable \
                and isinstance(observable[NOISE_DISTRIBUTION], str) \
                and observable[NOISE_DISTRIBUTION]:
            amici_val += observable[NOISE_DISTRIBUTION]
        else:
            amici_val += 'normal'
        amici_distrs[observable.name] = amici_val

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
    """
    Parse command line arguments

    :return:
        Parsed CLI arguments from :mod:`argparse`.
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
    parser.add_argument('-b', '--observables', dest='observable_file_name',
                        help='Observable table')

    parser.add_argument('-y', '--yaml', dest='yaml_file_name',
                       help='PEtab YAML problem filename')

    parser.add_argument('-n', '--model-name', dest='model_name',
                        help='Name of the python module generated for the '
                             'model')

    args = parser.parse_args()

    if not args.yaml_file_name \
            and not all((args.sbml_file_name, args.condition_file_name,
                         args.observable_file_name)):
        parser.error('When not specifying a model name or YAML file, then '
                     'SBML, condition and observable file must be specified')

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
            parameter_file=args.parameter_file_name,
            observable_files=args.observable_file_name)

    # First check for valid PEtab
    petab.lint_problem(pp)

    import_model(model_name=args.model_name,
                 sbml_model=pp.sbml_model,
                 condition_table=pp.condition_df,
                 observable_table=pp.observable_df,
                 model_output_dir=args.model_output_dir,
                 compile=args.compile,
                 verbose=args.verbose)


if __name__ == '__main__':
    main()
