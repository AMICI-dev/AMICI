"""
PySB-PEtab Import
------------
Import a model in the PySB-adapted :mod:`petab`
(https://github.com/PEtab-dev/PEtab) format into AMICI.
"""

import logging
import os
from typing import List, Dict, Union, Optional, Tuple, Iterable

import libsbml
import petab
import sympy as sp
from petab.C import (CONDITION_NAME, OBSERVABLE_TRANSFORMATION, LIN,
                     OBSERVABLE_FORMULA, NOISE_FORMULA, FORMAT_VERSION,
                     PARAMETER_FILE, SBML_FILES, CONDITION_FILES,
                     MEASUREMENT_FILES, VISUALIZATION_FILES, OBSERVABLE_FILES)

from .logging import get_logger, log_execution_time, set_log_level
from . import petab_import

logger = get_logger(__name__, logging.WARNING)


class PysbPetabProblem(petab.Problem):
    """Representation of a PySB-model-based PEtab problem

    This class extends :class:`petab.Problem` with a PySB model.
    The model is augmented with the observation model based on the PEtab
    observable table.
    For now, a dummy SBML model is created which allows used the existing
    SBML-PEtab API.

    :ivar pysb_model:
        PySB model instance from of this PEtab problem.

    """

    def __init__(self, pysb_model: 'pysb.Model' = None, *args, **kwargs):
        """
        Constructor

        :param pysb_model: PySB model instance for this PEtab problem
        :param args: See petab.Problem.__init__
        :param kwargs: See petab.Problem.__init__
        """

        super().__init__(*args, **kwargs)

        self.pysb_model: 'pysb.Model' = pysb_model
        self._add_observation_model()

        if self.pysb_model is not None:
            self.sbml_document, self.sbml_model = \
                create_dummy_sbml(
                    self.pysb_model,
                    observable_ids=self.observable_df.index.values
                    if self.observable_df is not None else None
                )

    def _add_observation_model(self):
        """Extend PySB model by observation model as defined in the PEtab
        observables table"""

        import pysb

        # add any required output parameters
        local_syms = {sp.Symbol.__str__(comp): comp for comp in
                      self.pysb_model.components if
                      isinstance(comp, sp.Symbol)}
        for formula in [*self.observable_df[OBSERVABLE_FORMULA],
                        *self.observable_df[NOISE_FORMULA]]:
            sym = sp.sympify(formula, locals=local_syms)
            for s in sym.free_symbols:
                if not isinstance(s, pysb.Component):
                    p = pysb.Parameter(str(s), 1.0)
                    local_syms[sp.Symbol.__str__(p)] = p

        # add observables and sigmas to pysb model
        for (observable_id, observable_formula, noise_formula) \
                in zip(self.observable_df.index,
                       self.observable_df[OBSERVABLE_FORMULA],
                       self.observable_df[NOISE_FORMULA]):
            # No observableTransformation so far
            if OBSERVABLE_TRANSFORMATION in self.observable_df:
                trafo = self.observable_df.loc[observable_id,
                                               OBSERVABLE_TRANSFORMATION]
                if trafo and trafo != LIN:
                    raise NotImplementedError(
                        "Observable transformation currently unsupported "
                        "for PySB models")
            obs_symbol = sp.sympify(observable_formula, locals=local_syms)
            obs_expr = pysb.Expression(observable_id, obs_symbol)
            local_syms[observable_id] = obs_expr

            sigma_id = f"{observable_id}_sigma"
            sigma_symbol = sp.sympify(
                noise_formula,
                locals=local_syms
            )
            sigma_expr = pysb.Expression(sigma_id, sigma_symbol)
            local_syms[sigma_id] = sigma_expr

    @staticmethod
    def from_files(condition_file: str = None,
                   measurement_file: Union[str, Iterable[str]] = None,
                   parameter_file: Union[str, List[str]] = None,
                   visualization_files: Union[str, Iterable[str]] = None,
                   observable_files: Union[str, Iterable[str]] = None,
                   pysb_model_file: str = None,
                   ) -> 'PysbPetabProblem':
        """
        Factory method to load model and tables from files.

        Arguments:
            condition_file: PEtab condition table
            measurement_file: PEtab measurement table
            parameter_file: PEtab parameter table
            visualization_files: PEtab visualization tables
            observable_files: PEtab observables tables
            pysb_model_file: PySB model file
        """

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
            pysb_model=pysb_model_from_path(
                pysb_model_file=pysb_model_file),
            condition_df=condition_df,
            measurement_df=measurement_df,
            parameter_df=parameter_df,
            observable_df=observable_df,
            visualization_df=visualization_df)

    @staticmethod
    def from_yaml(yaml_config: Union[Dict, str]) -> 'PysbPetabProblem':
        """
        Factory method to load model and tables as specified by YAML file.

        NOTE: The PySB model is currently expected in the YAML file under
        ``sbml_files``.

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


def create_dummy_sbml(
        pysb_model: 'pysb.Model',
        observable_ids: Optional[Iterable[str]] = None
) -> Tuple['libsbml.Model', 'libsbml.SBMLDocument']:
    """Create SBML dummy model for to use PySB models with PEtab.

    Model must at least contain PEtab problem parameter and noise parameters
    for observables.

    :param pysb_model: PySB model
    :param observable_ids: Observable IDs
    :return: A dummy SBML model and document.
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


@log_execution_time('Importing PEtab model', logger)
def import_model_pysb(
        petab_problem: PysbPetabProblem,
        model_output_dir: Optional[str] = None,
        verbose: Optional[Union[bool, int]] = True,
        **kwargs
) -> None:
    """
    Create AMICI model from PySB-PEtab problem

    :param petab_problem:
        PySB PEtab problem

    :param model_output_dir:
        Directory to write the model code to. Will be created if doesn't
        exist. Defaults to current directory.

    :param verbose:
        Print/log extra information.

    :param kwargs:
        Additional keyword arguments to be passed to
        :meth:`amici.pysb_import.pysb2amici`.
    """

    set_log_level(logger, verbose)

    logger.info(f"Importing model ...")

    observable_table = petab_problem.observable_df
    pysb_model = petab_problem.pysb_model

    # For pysb, we only allow parameters in the condition table
    # those must be pysb model parameters (either natively, or output
    # parameters from measurement or condition table that have been added in
    # PysbPetabProblem
    model_parameters = [p.name for p in pysb_model.parameters]
    for x in petab_problem.condition_df.columns:
        if x == CONDITION_NAME:
            continue

        if x not in model_parameters:
            raise NotImplementedError(
                "For PySB PEtab import, only model parameters, but no states "
                "or compartments are allowed in the condition table."
                f"Offending column: {x}"
            )

    constant_parameters = petab_import.get_fixed_parameters(
        petab_problem.sbml_model, petab_problem.condition_df)

    if observable_table is None:
        observables = None
        sigmas = None
    else:
        observables = [expr.name for expr in pysb_model.expressions
                       if expr.name in observable_table.index]

        sigmas = {obs_id: f"{obs_id}_sigma" for obs_id in observables}

    from amici.pysb_import import pysb2amici
    pysb2amici(pysb_model, model_output_dir, verbose=True,
               observables=observables,
               sigmas=sigmas,
               constant_parameters=constant_parameters,
               **kwargs)
