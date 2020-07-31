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
from _collections import OrderedDict
from itertools import chain
from typing import List, Dict, Union, Optional, Tuple, Iterable

import amici
import libsbml
import pandas as pd
import petab
import sympy as sp

from amici.logging import get_logger, log_execution_time, set_log_level
from petab.C import *
from amici.pysb_import import pysb2amici

logger = get_logger(__name__, logging.WARNING)

# ID of model parameter that is to be added to SBML model to indicate
#  preequilibration
PREEQ_INDICATOR_ID = 'preequilibration_indicator'


class PysbPetabProblem(petab.Problem):
    def __init__(self, pysb_model_module: 'pysb.Model' = None,
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.pysb_model_module = pysb_model_module
        import pysb

        # add any required output parameters
        # sp.Symbol.__str__(comp)
        locals = {sp.Symbol.__str__(comp): comp for comp in
                  pysb_model_module.components if
                  isinstance(comp, sp.Symbol)}
        for formula in [*self.observable_df[petab.OBSERVABLE_FORMULA],
                       *self.observable_df[petab.NOISE_FORMULA]]:
            sym = sp.sympify(formula, locals=locals)
            for s in sym.free_symbols:
                if not isinstance(s, pysb.Component):
                    p = pysb.Parameter(str(s), 1.0)
                    locals[sp.Symbol.__str__(p)] = p


        # add observables and sigmas to pysb model
        for (observable_id, observable_formula, noise_formula) \
                in zip(self.observable_df.index,
                       self.observable_df[petab.OBSERVABLE_FORMULA],
                       self.observable_df[petab.NOISE_FORMULA]):
            if petab.OBSERVABLE_TRANSFORMATION in self.observable_df:
                trafo = self.observable_df.loc[observable_id, petab.OBSERVABLE_TRANSFORMATION]
                if trafo and trafo != petab.LIN:
                    raise NotImplementedError("Observable transformation currently unsupported for PySB models")
            obs_symbol = sp.sympify(
                    observable_formula,
                    locals=locals
                )
            ex = pysb.Expression(observable_id, obs_symbol)
            locals[observable_id] = ex

            sigma_id = f"{observable_id}_sigma"  # TODO naming?
            sigma_symbol = sp.sympify(
                    noise_formula,
                    locals=locals
                )
            ex = pysb.Expression(sigma_id, sigma_symbol)
            locals[sigma_id] = ex

        if self.pysb_model_module is not None:
            self.sbml_document, self.sbml_model = create_dummy_sbml(
                self.pysb_model_module,
                observable_ids=self.observable_df.index.values
                if self.observable_df is not None else None
            )

    @staticmethod
    def from_files(sbml_file: str = None,
                   condition_file: str = None,
                   measurement_file: Union[str, Iterable[str]] = None,
                   parameter_file: Union[str, List[str]] = None,
                   visualization_files: Union[str, Iterable[str]] = None,
                   observable_files: Union[str, Iterable[str]] = None,
                   pysb_model_file: str = None,
                   ) -> 'PysbPetabProblem':
        """
        Factory method to load model and tables from files.

        Arguments:
            sbml_file: PEtab SBML model
            condition_file: PEtab condition table
            measurement_file: PEtab measurement table
            parameter_file: PEtab parameter table
            visualization_files: PEtab visualization tables
            observable_files: PEtab observables tables
            pysb_model_file: PySB model file
        """

        sbml_model = sbml_document = sbml_reader = None
        condition_df = measurement_df = parameter_df = visualization_df = None
        observable_df = None

        if condition_file:
            condition_df = petab.conditions.get_condition_df(condition_file)

        if measurement_file:
            # If there are multiple tables, we will merge them
            measurement_df = petab.core.concat_tables(
                measurement_file, petab.measurements.get_measurement_df)

        if parameter_file:
            parameter_df = petab.parameters.get_parameter_df(parameter_file)

        if sbml_file:
            sbml_reader = libsbml.SBMLReader()
            sbml_document = sbml_reader.readSBML(sbml_file)
            sbml_model = sbml_document.getModel()

        if visualization_files:
            # If there are multiple tables, we will merge them
            visualization_df = petab.core.concat_tables(
                visualization_files, petab.core.get_visualization_df)

        if observable_files:
            # If there are multiple tables, we will merge them
            observable_df = petab.core.concat_tables(
                observable_files, petab.observables.get_observable_df)
        from amici.pysb_import import pysb_model_from_path
        return PysbPetabProblem(
            pysb_model_module=pysb_model_from_path(
                pysb_model_file=pysb_model_file),
            condition_df=condition_df,
            measurement_df=measurement_df,
            parameter_df=parameter_df,
            observable_df=observable_df,
            sbml_model=sbml_model,
            sbml_document=sbml_document,
            sbml_reader=sbml_reader,
            visualization_df=visualization_df)

    @staticmethod
    def from_yaml(yaml_config: Union[Dict, str]) -> 'PysbPetabProblem':
        """
        Factory method to load model and tables as specified by YAML file.

        Arguments:
            yaml_config: PEtab configuration as dictionary or YAML file name
        """
        from petab.yaml import (load_yaml, is_composite_problem,
                                assert_single_condition_and_sbml_file)
        if isinstance(yaml_config, str):
            path_prefix = os.path.dirname(yaml_config)
            yaml_config = load_yaml(yaml_config)
        else:
            path_prefix = ""

        if is_composite_problem(yaml_config):
            raise ValueError('petab.Problem.from_yaml() can only be used for '
                             'yaml files comprising a single model. '
                             'Consider using '
                             'petab.CompositeProblem.from_yaml() instead.')

        if yaml_config[FORMAT_VERSION] != petab.__format_version__:
            raise ValueError("Provided PEtab files are of unsupported version"
                             f"{yaml_config[FORMAT_VERSION]}. Expected "
                             f"{petab.__format_version__}.")

        problem0 = yaml_config['problems'][0]

        assert_single_condition_and_sbml_file(problem0)

        if isinstance(yaml_config[PARAMETER_FILE], list):
            parameter_file = [
                os.path.join(path_prefix, f)
                for f in yaml_config[PARAMETER_FILE]
            ]
        else:
            parameter_file = os.path.join(
                path_prefix, yaml_config[PARAMETER_FILE])

        return PysbPetabProblem.from_files(
            pysb_model_file=os.path.join(
                path_prefix, problem0[SBML_FILES][0]),
            measurement_file=[os.path.join(path_prefix, f)
                              for f in problem0[MEASUREMENT_FILES]],
            condition_file=os.path.join(
                path_prefix, problem0[CONDITION_FILES][0]),
            parameter_file=parameter_file,
            visualization_files=[
                os.path.join(path_prefix, f)
                for f in problem0.get(VISUALIZATION_FILES, [])],
            observable_files=[
                os.path.join(path_prefix, f)
                for f in problem0.get(OBSERVABLE_FILES, [])]
        )


def create_dummy_sbml(pysb_model: 'pysb.Model', observable_ids=None
                      ) -> Tuple['libsbml.Model', 'libsbml.SBMLDocument']:
    """Create SBML dummy model for to use pysb Models with PEtab.

    Model must at least contain petab problem parameter and noise parameters
    for observables.
    """

    import libsbml

    document = libsbml.SBMLDocument(3, 1)
    dummy_sbml_model = document.createModel()
    dummy_sbml_model.setTimeUnits("second")
    dummy_sbml_model.setExtentUnits("mole")
    dummy_sbml_model.setSubstanceUnits('mole')

    for species in pysb_model.parameters:
        p = dummy_sbml_model.createParameter()
        p.setId(species.name)
        p.setConstant(True)
        p.setValue(0.0)

    for observable_id in observable_ids:
        p = dummy_sbml_model.createParameter()
        p.setId(f"noiseParameter1_{observable_id}")
        p.setConstant(True)
        p.setValue(0.0)

    return document, dummy_sbml_model


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
                # remove from reactants
                pass
            while reaction.removeProduct(species_id):
                # remove from products
                pass
            while reaction.removeModifier(species_id):
                # remove from modifiers
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
        if isinstance(petab_problem, PysbPetabProblem):
            raise ValueError("Parameter `model_output_dir` is required.")

        model_output_dir = \
            _create_model_output_dir_name(petab_problem.sbml_model)
    else:
        model_output_dir = os.path.abspath(model_output_dir)

    if isinstance(petab_problem, PysbPetabProblem):
        if model_name is None:
            model_name = petab_problem.pysb_model_module.name
        else:
            raise ValueError(
                "Argument model_name currently not allowed for pysb models")
    elif model_name is None:
        model_name = _create_model_name(model_output_dir)

    # create folder
    if not os.path.exists(model_output_dir):
        os.makedirs(model_output_dir)

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
        if isinstance(petab_problem, PysbPetabProblem):
            import_model_pysb(
                petab_problem,
                model_output_dir=model_output_dir,
                **kwargs)
        else:
            import_model(
                sbml_model=petab_problem.sbml_model,
                condition_table=petab_problem.condition_df,
                observable_table=petab_problem.observable_df,
                measurement_table=petab_problem.measurement_df,
                model_name=model_name,
                model_output_dir=model_output_dir,
                **kwargs)

    # import model
    model_module = amici.import_model_module(model_name, model_output_dir)
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
        model_module = importlib.import_module(model_name)
    except ModuleNotFoundError:
        return False

    # no need to (re-)compile
    return hasattr(model_module, "getModel")


@log_execution_time('Importing PEtab model', logger)
def import_model(sbml_model: Union[str, 'libsbml.Model'],
                 condition_table: Optional[Union[str, pd.DataFrame]] = None,
                 observable_table: Optional[Union[str, pd.DataFrame]] = None,
                 measurement_table: Optional[Union[str, pd.DataFrame]] = None,
                 model_name: Optional[str] = None,
                 model_output_dir: Optional[str] = None,
                 verbose: Optional[Union[bool, int]] = True,
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

    :param measurement_table:
        PEtab measurement table.

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
    formulas = chain((val['formula'] for val in observables.values()),
                     sigmas.values())
    output_parameters = OrderedDict()
    for formula in formulas:
        # we want reproducible parameter ordering upon repeated import
        free_syms = sorted(sp.sympify(formula).free_symbols,
                           key=lambda symbol: symbol.name)
        for free_sym in free_syms:
            sym = str(free_sym)
            if sbml_model.getElementBySId(sym) is None and sym != 'time':
                output_parameters[sym] = None
    logger.debug(f"Adding output parameters to model: {output_parameters}")
    for par in output_parameters.keys():
        petab.add_global_parameter(sbml_model, par)
    # <EndWorkAround>

    # TODO: to parameterize initial states or compartment sizes, we currently
    #  need initial assignments. if they occur in the condition table, we
    #  create a new parameter initial_${startOrCompartmentID}.
    #  feels dirty and should be changed (see also #924)
    # <BeginWorkAround>
    initial_states = [col for col in condition_df
                      if sbml_model.getSpecies(col) is not None]
    initial_sizes = [col for col in condition_df
                     if sbml_model.getCompartment(col) is not None]
    fixed_parameters = []
    if len(initial_states) or len(initial_sizes):
        # add preequilibration indicator variable
        # NOTE: would only be required if we actually have preequilibration
        #  adding it anyways. can be optimized-out later
        if sbml_model.getParameter(PREEQ_INDICATOR_ID) is not None:
            raise AssertionError("Model already has a parameter with ID "
                                 f"{PREEQ_INDICATOR_ID}. Cannot handle "
                                 "species and compartments in condition table "
                                 "then.")
        indicator = sbml_model.createParameter()
        indicator.setId(PREEQ_INDICATOR_ID)
        indicator.setName(PREEQ_INDICATOR_ID)
        # Can only reset parameters after preequilibration if they are fixed.
        fixed_parameters.append(PREEQ_INDICATOR_ID)

    for assignee_id in initial_sizes + initial_states:
        init_par_id_preeq = f"initial_{assignee_id}_preeq"
        init_par_id_sim = f"initial_{assignee_id}_sim"
        for init_par_id in [init_par_id_preeq, init_par_id_sim]:
            if sbml_model.getElementBySId(init_par_id) is not None:
                raise ValueError(
                    "Cannot create parameter for initial assignment "
                    f"for {assignee_id} because an entity named "
                    f"{init_par_id} exists already in the model.")
            init_par = sbml_model.createParameter()
            init_par.setId(init_par_id)
            init_par.setName(init_par_id)
        assignment = sbml_model.createInitialAssignment()
        assignment.setSymbol(assignee_id)
        formula = f'{PREEQ_INDICATOR_ID} * {init_par_id_preeq} ' \
                  f'+ (1 - {PREEQ_INDICATOR_ID}) * {init_par_id_sim}'
        math_ast = libsbml.parseL3Formula(formula)
        assignment.setMath(math_ast)
    # <EndWorkAround>

    fixed_parameters.extend(
        get_fixed_parameters(sbml_model=sbml_model, condition_df=condition_df))

    logger.debug(f"Fixed parameters are {fixed_parameters}")
    logger.info(f"Overall fixed parameters: {len(fixed_parameters)}")
    logger.info("Variable parameters: "
                + str(len(sbml_model.getListOfParameters())
                      - len(fixed_parameters)))

    # Create Python module from SBML model
    sbml_importer.sbml2amici(
        model_name=model_name,
        output_dir=model_output_dir,
        observables=observables,
        constant_parameters=fixed_parameters,
        sigmas=sigmas,
        allow_reinit_fixpar_initcond=allow_reinit_fixpar_initcond,
        noise_distributions=noise_distrs,
        verbose=verbose,
        **kwargs)


import_model_sbml = import_model


@log_execution_time('Importing PEtab model', logger)
def import_model_pysb(
        petab_problem: PysbPetabProblem,
        model_output_dir: Optional[str] = None,
        verbose: Optional[Union[bool, int]] = True,
        **kwargs
) -> None:
    """
    Create AMICI model from PEtab problem

    :param petab_problem:
        PySB PEtab problem

    :param model_output_dir:
        Directory to write the model code to. Will be created if doesn't
        exist. Defaults to current directory.

    :param verbose:
        Print/log extra information.

    :param kwargs:
        Additional keyword arguments to be passed to
        :meth:`amici.sbml_import.SbmlImporter.sbml2amici`.
    """

    set_log_level(logger, verbose)

    logger.info(f"Importing model ...")

    observable_table = petab_problem.observable_df
    pysb_model_module = petab_problem.pysb_model_module

    # For pysb, we only allow parameters in the condition table
    # those must be pysb model parameters (either natively, or output
    # parameters from measurement or condition table that have been added in
    # PysbPetabProblem
    model_parameters = [p.name for p in pysb_model_module.parameters]
    for x in petab_problem.condition_df.columns:
        if x == petab.CONDITION_NAME: continue
        if x not in model_parameters:
            raise NotImplementedError(
                "For PySB PEtab import, only model parameters, but no states "
                "or compartments are allowed in the condition table."
                f"Offending column: {x}"
            )


    constant_parameters = get_fixed_parameters(petab_problem.sbml_model,
                                               petab_problem.condition_df)

    if observable_table is None:
        observables = None
        sigmas = None
    else:
        observables = [expr.name for expr in pysb_model_module.expressions
                       if expr.name in observable_table.index]

        def get_expr(x):
            for expr in pysb_model_module.expressions:
                if expr.name == x:
                    return expr
        sigmas = {obs_id: get_expr(f"{obs_id}_sigma") for obs_id in observables}

    pysb2amici(pysb_model_module, model_output_dir, verbose=True,
               observables=observables,
               sigmas=sigmas,
               constant_parameters=constant_parameters,
               # compute_conservation_laws=False,
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
                 measurement_table=pp.measurement_df,
                 model_output_dir=args.model_output_dir,
                 compile=args.compile,
                 verbose=args.verbose)


if __name__ == '__main__':
    main()
