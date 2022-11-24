"""
PEtab Import
------------
Import a model in the :mod:`petab` (https://github.com/PEtab-dev/PEtab) format
into AMICI.
"""
import argparse
import importlib
import logging
import os
import re
import shutil
import tempfile
from _collections import OrderedDict
from itertools import chain
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from warnings import warn

import libsbml
import pandas as pd
import petab
import sympy as sp
from petab.C import *
from petab.parameters import get_valid_parameters_for_parameter_table

import amici
from amici.logging import get_logger, log_execution_time, set_log_level

try:
    from amici.petab_import_pysb import PysbPetabProblem, import_model_pysb
except ModuleNotFoundError:
    # pysb not available
    PysbPetabProblem = None
    import_model_pysb = None

logger = get_logger(__name__, logging.WARNING)

# ID of model parameter that is to be added to SBML model to indicate
#  preequilibration
PREEQ_INDICATOR_ID = 'preequilibration_indicator'


def _add_global_parameter(sbml_model: libsbml.Model,
                          parameter_id: str,
                          parameter_name: str = None,
                          constant: bool = False,
                          units: str = 'dimensionless',
                          value: float = 0.0) -> libsbml.Parameter:
    """Add new global parameter to SBML model

    Arguments:
        sbml_model: SBML model
        parameter_id: ID of the new parameter
        parameter_name: Name of the new parameter
        constant: Is parameter constant?
        units: SBML unit ID
        value: parameter value

    Returns:
        The created parameter
    """
    if parameter_name is None:
        parameter_name = parameter_id

    p = sbml_model.createParameter()
    p.setId(parameter_id)
    p.setName(parameter_name)
    p.setConstant(constant)
    p.setValue(value)
    p.setUnits(units)
    return p


def get_fixed_parameters(
        petab_problem: petab.Problem,
        non_estimated_parameters_as_constants=True,
) -> List[str]:
    """
    Determine, set and return fixed model parameters.

    Non-estimated parameters and parameters specified in the condition table
    are turned into constants (unless they are overridden).
    Only global SBML parameters are considered. Local parameters are ignored.

    :param petab_problem:
        The PEtab problem instance

    :param non_estimated_parameters_as_constants:
        Whether parameters marked as non-estimated in PEtab should be
        considered constant in AMICI. Setting this to ``True`` will reduce
        model size and simulation times. If sensitivities with respect to those
        parameters are required, this should be set to ``False``.

    :return:
        List of IDs of parameters which are to be considered constant.
    """
    # initial concentrations for species or initial compartment sizes in
    # condition table will need to be turned into fixed parameters

    # if there is no initial assignment for that species, we'd need
    # to create one. to avoid any naming collision right away, we don't
    # allow that for now

    # we can't handle them yet
    compartments = [
        col for col in petab_problem.condition_df
        if petab_problem.sbml_model.getCompartment(col) is not None
    ]
    if compartments:
        raise NotImplementedError("Can't handle initial compartment sizes "
                                  "at the moment. Consider creating an "
                                  f"initial assignment for {compartments}")

    # if we have a parameter table, all parameters that are allowed to be
    #  listed in the parameter table, but are not marked as estimated, can be
    #  turned into AMICI constants
    # due to legacy API, we might not always have a parameter table, though
    fixed_parameters = set()
    if petab_problem.parameter_df is not None:
        all_parameters = get_valid_parameters_for_parameter_table(
            model=petab_problem.model,
            condition_df=petab_problem.condition_df,
            observable_df=petab_problem.observable_df
            if petab_problem.observable_df is not None
            else pd.DataFrame(columns=petab.OBSERVABLE_DF_REQUIRED_COLS),
            measurement_df=petab_problem.measurement_df
            if petab_problem.measurement_df is not None
            else pd.DataFrame(columns=petab.MEASUREMENT_DF_REQUIRED_COLS),
        )
        if non_estimated_parameters_as_constants:
            estimated_parameters = \
                petab_problem.parameter_df.index.values[
                    petab_problem.parameter_df[ESTIMATE] == 1]
        else:
            # don't treat parameter table parameters as constants
            estimated_parameters = petab_problem.parameter_df.index.values
        fixed_parameters = set(all_parameters) - set(estimated_parameters)

    sbml_model = petab_problem.sbml_model
    condition_df = petab_problem.condition_df

    # Column names are model parameter IDs, compartment IDs or species IDs.
    # Thereof, all parameters except for any overridden ones should be made
    # constant.
    # (Could potentially still be made constant, but leaving them might
    # increase model reusability)

    # handle parameters in condition table
    if condition_df is not None:
        logger.debug(f'Condition table: {condition_df.shape}')

        # remove overridden parameters (`object`-type columns)
        fixed_parameters.update(
            p for p in condition_df.columns
            # get rid of conditionName column
            if p != CONDITION_NAME
            # there is no parametric override
            # TODO: could check if the final overriding parameter is estimated
            #  or not, but for now, we skip the parameter if there is any kind
            #  of overriding
            if condition_df[p].dtype != 'O'
               # p is a parameter
               and sbml_model.getParameter(p) is not None
               # but not a rule target
               and sbml_model.getRuleByVariable(p) is None
        )

    # Ensure mentioned parameters exist in the model. Remove additional ones
    # from list
    for fixed_parameter in fixed_parameters.copy():
        # check global parameters
        if not sbml_model.getParameter(fixed_parameter):
            logger.warning(f"Parameter or species '{fixed_parameter}'"
                           " provided in condition table but not present in"
                           " model. Ignoring.")
            fixed_parameters.remove(fixed_parameter)

    return list(sorted(fixed_parameters))


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


def import_petab_problem(
        petab_problem: petab.Problem,
        model_output_dir: Union[str, Path, None] = None,
        model_name: str = None,
        force_compile: bool = False,
        non_estimated_parameters_as_constants = True,
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
        the model ID will be used.

    :param force_compile:
        Whether to compile the model even if the target folder is not empty,
        or the model exists already.

    :param non_estimated_parameters_as_constants:
        Whether parameters marked as non-estimated in PEtab should be
        considered constant in AMICI. Setting this to ``True`` will reduce
        model size and simulation times. If sensitivities with respect to those
        parameters are required, this should be set to ``False``.

    :param kwargs:
        Additional keyword arguments to be passed to
        :meth:`amici.sbml_import.SbmlImporter.sbml2amici`.

    :return:
        The imported model.
    """
    # generate folder and model name if necessary
    if model_output_dir is None:
        if PysbPetabProblem and isinstance(petab_problem, PysbPetabProblem):
            raise ValueError("Parameter `model_output_dir` is required.")

        model_output_dir = \
            _create_model_output_dir_name(petab_problem.sbml_model)
    else:
        model_output_dir = os.path.abspath(model_output_dir)

    if PysbPetabProblem and isinstance(petab_problem, PysbPetabProblem) \
            and model_name is None:
        model_name = petab_problem.pysb_model.name
    elif model_name is None:
        model_name = _create_model_name(model_output_dir)

    # create folder
    if not os.path.exists(model_output_dir):
        os.makedirs(model_output_dir)

    # check if compilation necessary
    if force_compile or not _can_import_model(model_name, model_output_dir):
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
        if PysbPetabProblem and isinstance(petab_problem, PysbPetabProblem):
            import_model_pysb(
                petab_problem,
                model_name=model_name,
                model_output_dir=model_output_dir,
                **kwargs)
        else:
            import_model_sbml(
                petab_problem=petab_problem,
                model_name=model_name,
                model_output_dir=model_output_dir,
                non_estimated_parameters_as_constants=
                non_estimated_parameters_as_constants,
                **kwargs)

    # import model
    model_module = amici.import_model_module(model_name, model_output_dir)
    model = model_module.getModel()
    check_model(amici_model=model, petab_problem=petab_problem)

    logger.info(f"Successfully loaded model {model_name} "
                f"from {model_output_dir}.")

    return model


def check_model(
    amici_model: amici.Model,
    petab_problem: petab.Problem,
) -> None:
    """Check that the model is consistent with the PEtab problem."""
    if petab_problem.parameter_df is None:
        return

    amici_ids_free = set(amici_model.getParameterIds())
    amici_ids = amici_ids_free | set(amici_model.getFixedParameterIds())

    petab_ids_free = set(petab_problem.parameter_df.loc[
        petab_problem.parameter_df[ESTIMATE] == 1
    ].index)

    amici_ids_free_required = petab_ids_free.intersection(amici_ids)

    if not amici_ids_free_required.issubset(amici_ids_free):
        raise ValueError(
            'The available AMICI model does not support estimating the '
            'following parameters. Please recompile the model and ensure '
            'that these parameters are not treated as constants. Deleting '
            'the current model might also resolve this. Parameters: '
            f'{amici_ids_free_required.difference(amici_ids_free)}'
        )


def _create_model_output_dir_name(sbml_model: 'libsbml.Model') -> Path:
    """
    Find a folder for storing the compiled amici model.
    If possible, use the sbml model id, otherwise create a random folder.
    The folder will be located in the `amici_models` subfolder of the current
    folder.
    """
    BASE_DIR = Path("amici_models").absolute()
    BASE_DIR.mkdir(exist_ok=True)
    # try sbml model id
    if sbml_model_id := sbml_model.getId():
        return BASE_DIR / sbml_model_id

    # create random folder name
    return Path(tempfile.mkdtemp(dir=BASE_DIR))


def _create_model_name(folder: Union[str, Path]) -> str:
    """
    Create a name for the model.
    Just re-use the last part of the folder.
    """
    return os.path.split(os.path.normpath(folder))[-1]


def _can_import_model(
        model_name: str,
        model_output_dir: Union[str, Path]
) -> bool:
    """
    Check whether a module of that name can already be imported.
    """
    # try to import (in particular checks version)
    try:
        with amici.add_path(model_output_dir):
            model_module = importlib.import_module(model_name)
    except ModuleNotFoundError:
        return False

    # no need to (re-)compile
    return hasattr(model_module, "getModel")


@log_execution_time('Importing PEtab model', logger)
def import_model_sbml(
        sbml_model: Union[str, Path, 'libsbml.Model'] = None,
        condition_table: Optional[Union[str, Path, pd.DataFrame]] = None,
        observable_table: Optional[Union[str, Path, pd.DataFrame]] = None,
        measurement_table: Optional[Union[str, Path, pd.DataFrame]] = None,
        petab_problem: petab.Problem = None,
        model_name: Optional[str] = None,
        model_output_dir: Optional[Union[str, Path]] = None,
        verbose: Optional[Union[bool, int]] = True,
        allow_reinit_fixpar_initcond: bool = True,
        validate: bool = True,
        non_estimated_parameters_as_constants=True,
        **kwargs) -> amici.SbmlImporter:
    """
    Create AMICI model from PEtab problem

    :param sbml_model:
        PEtab SBML model or SBML file name.
        Deprecated, pass ``petab_problem`` instead.

    :param condition_table:
        PEtab condition table. If provided, parameters from there will be
        turned into AMICI constant parameters (i.e. parameters w.r.t. which
        no sensitivities will be computed).
        Deprecated, pass ``petab_problem`` instead.

    :param observable_table:
        PEtab observable table. Deprecated, pass ``petab_problem`` instead.

    :param measurement_table:
        PEtab measurement table. Deprecated, pass ``petab_problem`` instead.

    :param petab_problem:
        PEtab problem.

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

    :param validate:
        Whether to validate the PEtab problem

    :param non_estimated_parameters_as_constants:
        Whether parameters marked as non-estimated in PEtab should be
        considered constant in AMICI. Setting this to ``True`` will reduce
        model size and simulation times. If sensitivities with respect to those
        parameters are required, this should be set to ``False``.

    :param kwargs:
        Additional keyword arguments to be passed to
        :meth:`amici.sbml_import.SbmlImporter.sbml2amici`.

    :return:
        The created :class:`amici.sbml_import.SbmlImporter` instance.
    """
    from petab.models.sbml_model import SbmlModel

    set_log_level(logger, verbose)

    logger.info("Importing model ...")

    if any([sbml_model, condition_table, observable_table, measurement_table]):
        warn("The `sbml_model`, `condition_table`, `observable_table`, and "
             "`measurement_table` arguments are deprecated and will be "
             "removed in a future version. Use `petab_problem` instead.",
             DeprecationWarning, stacklevel=2)
        if petab_problem:
            raise ValueError("Must not pass a `petab_problem` argument in "
                             "combination with any of `sbml_model`, "
                             "`condition_table`, `observable_table`, or "
                             "`measurement_table`.")

        petab_problem = petab.Problem(
            model=SbmlModel(sbml_model)
            if isinstance(sbml_model, libsbml.Model)
            else SbmlModel.from_file(sbml_model),
            condition_df=petab.get_condition_df(condition_table),
            observable_df=petab.get_observable_df(observable_table),
        )

    if petab_problem.observable_df is None:
        raise NotImplementedError("PEtab import without observables table "
                                  "is currently not supported.")

    assert isinstance(petab_problem.model, SbmlModel)

    if validate:
        logger.info("Validating PEtab problem ...")
        petab.lint_problem(petab_problem)

    # Model name from SBML ID or filename
    if model_name is None:
        if not (model_name := petab_problem.model.sbml_model.getId()):
            if not isinstance(sbml_model, (str, Path)):
                raise ValueError("No `model_name` was provided and no model "
                                 "ID was specified in the SBML model.")
            model_name = os.path.splitext(os.path.split(sbml_model)[-1])[0]

    if model_output_dir is None:
        model_output_dir = os.path.join(
            os.getcwd(), f"{model_name}-amici{amici.__version__}"
        )

    logger.info(f"Model name is '{model_name}'.\n"
                f"Writing model code to '{model_output_dir}'.")

    # Create a copy, because it will be modified by SbmlImporter
    sbml_doc = petab_problem.model.sbml_model.getSBMLDocument().clone()
    sbml_model = sbml_doc.getModel()

    show_model_info(sbml_model)

    sbml_importer = amici.SbmlImporter(sbml_model)
    sbml_model = sbml_importer.sbml

    allow_n_noise_pars = \
        not petab.lint.observable_table_has_nontrivial_noise_formula(
            petab_problem.observable_df
        )
    if petab_problem.measurement_df is not None and \
            petab.lint.measurement_table_has_timepoint_specific_mappings(
                petab_problem.measurement_df,
                allow_scalar_numeric_noise_parameters=allow_n_noise_pars
            ):
        raise ValueError(
            'AMICI does not support importing models with timepoint specific '
            'mappings for noise or observable parameters. Please flatten '
            'the problem and try again.'
        )

    if petab_problem.observable_df is not None:
        observables, noise_distrs, sigmas = \
            get_observation_model(petab_problem.observable_df)
    else:
        observables = noise_distrs = sigmas = None

    logger.info(f'Observables: {len(observables)}')
    logger.info(f'Sigmas: {len(sigmas)}')

    if len(sigmas) != len(observables):
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
            if sbml_model.getElementBySId(sym) is None and sym != 'time' \
                    and sym not in observables:
                output_parameters[sym] = None
    logger.debug("Adding output parameters to model: "
                 f"{list(output_parameters.keys())}")
    for par in output_parameters.keys():
        _add_global_parameter(sbml_model, par)
    # <EndWorkAround>

    # TODO: to parameterize initial states or compartment sizes, we currently
    #  need initial assignments. if they occur in the condition table, we
    #  create a new parameter initial_${startOrCompartmentID}.
    #  feels dirty and should be changed (see also #924)
    # <BeginWorkAround>

    initial_states = [col for col in petab_problem.condition_df
                      if element_is_state(sbml_model, col)]
    fixed_parameters = []
    if initial_states:
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
        logger.debug("Adding preequilibration indicator "
                     f"constant {PREEQ_INDICATOR_ID}")
    logger.debug(f"Adding initial assignments for {initial_states}")
    for assignee_id in initial_states:
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
        assignment = sbml_model.getInitialAssignment(assignee_id)
        if assignment is None:
            assignment = sbml_model.createInitialAssignment()
            assignment.setSymbol(assignee_id)
        else:
            logger.debug('The SBML model has an initial assignment defined '
                         f'for model entity {assignee_id}, but this entity '
                         'also has an initial value defined in the PEtab '
                         'condition table. The SBML initial assignment will '
                         'be overwritten to handle preequilibration and '
                         'initial values specified by the PEtab problem.')
        formula = f'{PREEQ_INDICATOR_ID} * {init_par_id_preeq} ' \
                  f'+ (1 - {PREEQ_INDICATOR_ID}) * {init_par_id_sim}'
        math_ast = libsbml.parseL3Formula(formula)
        assignment.setMath(math_ast)
    # <EndWorkAround>

    fixed_parameters.extend(
        get_fixed_parameters(
            petab_problem=petab_problem,
            non_estimated_parameters_as_constants=
            non_estimated_parameters_as_constants,
        ))

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

    if kwargs.get('compile', amici._get_default_argument(
            sbml_importer.sbml2amici, 'compile')):
        # check that the model extension was compiled successfully
        model_module = amici.import_model_module(model_name, model_output_dir)
        model = model_module.getModel()
        check_model(amici_model=model, petab_problem=petab_problem)

    return sbml_importer


# for backwards compatibility
import_model = import_model_sbml


def get_observation_model(
        observable_df: pd.DataFrame,
) -> Tuple[Dict[str, Dict[str, str]], Dict[str, str],
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

    nan_pat = r'^[nN]a[nN]$'
    for _, observable in observable_df.iterrows():
        oid = str(observable.name)
        # need to sanitize due to https://github.com/PEtab-dev/PEtab/issues/447
        name = re.sub(nan_pat, '', str(observable.get(OBSERVABLE_NAME, '')))
        formula_obs = re.sub(nan_pat, '', str(observable[OBSERVABLE_FORMULA]))
        formula_noise = re.sub(nan_pat, '', str(observable[NOISE_FORMULA]))
        observables[oid] = {'name': name, 'formula': formula_obs}
        sigmas[oid] = formula_noise

    # PEtab does currently not allow observables in noiseFormula and AMICI
    #  cannot handle states in sigma expressions. Therefore, where possible,
    #  replace species occurring in error model definition by observableIds.
    replacements = {
        sp.sympify(observable['formula']): sp.Symbol(observable_id)
        for observable_id, observable in observables.items()
    }
    for observable_id, formula in sigmas.items():
        repl = sp.sympify(formula).subs(replacements)
        sigmas[observable_id] = str(repl)

    noise_distrs = petab_noise_distributions_to_amici(observable_df)

    return observables, noise_distrs, sigmas


def petab_noise_distributions_to_amici(observable_df: pd.DataFrame
                                       ) -> Dict[str, str]:
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


def element_is_state(sbml_model: libsbml.Model, sbml_id: str) -> bool:
    """Does the element with ID `sbml_id` correspond to a state variable?
    """
    if sbml_model.getCompartment(sbml_id) is not None:
        return True
    if sbml_model.getSpecies(sbml_id) is not None:
        return True
    if (rule := sbml_model.getRuleByVariable(sbml_id)) is not None \
            and rule.getTypeCode() == libsbml.SBML_RATE_RULE:
        return True

    return False


def _parse_cli_args():
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
    parser.add_argument('--flatten', dest='flatten', default=False,
                        action='store_true',
                        help='Flatten measurement specific overrides of '
                             'observable and noise parameters')
    parser.add_argument('--no-sensitivities', dest='generate_sensitivity_code',
                        default=True, action='store_false',
                        help='Skip generation of sensitivity code')

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
    (https://github.com/PEtab-dev/PEtab/) format into AMICI.
    """
    args = _parse_cli_args()

    if args.yaml_file_name:
        pp = petab.Problem.from_yaml(args.yaml_file_name)
    else:
        pp = petab.Problem.from_files(
            sbml_file=args.sbml_file_name,
            condition_file=args.condition_file_name,
            measurement_file=args.measurement_file_name,
            parameter_file=args.parameter_file_name,
            observable_files=args.observable_file_name)

    # Check for valid PEtab before potentially modifying it
    petab.lint_problem(pp)

    if args.flatten:
        petab.flatten_timepoint_specific_output_overrides(pp)

    import_model(model_name=args.model_name,
                 sbml_model=pp.sbml_model,
                 condition_table=pp.condition_df,
                 observable_table=pp.observable_df,
                 measurement_table=pp.measurement_df,
                 model_output_dir=args.model_output_dir,
                 compile=args.compile,
                 generate_sensitivity_code=args.generate_sensitivity_code,
                 verbose=args.verbose,
                 validate=False)


if __name__ == '__main__':
    main()
