"""
Adapters for using AMICI with the `fiddy <https://github.com/ICB-DCM/fiddy/>`__
package for finite difference checks.


.. note::

    Like fiddy, this module is experimental and subject to change.
"""

from __future__ import annotations

import warnings
from collections.abc import Callable
from functools import partial
from inspect import signature
from typing import TYPE_CHECKING, Any

import numpy as np
import petab.v1 as petab
from fiddy import CachedFunction, Type, fiddy_array
from fiddy.directional_derivative import DirectionalDerivative
from fiddy.success import Consistency
from petab.v1.C import LIN, LOG, LOG10

from amici.sim.sundials import (
    AmiciExpData,
    AmiciModel,
    AmiciSolver,
    Model,
    ReturnData,
    SensitivityOrder,
    run_simulation,
)
from amici.sim.sundials.petab.v1 import LLH, SLLH, create_edatas
from amici.sim.sundials.petab.v1._parameter_mapping import (
    create_parameter_mapping,
)

if TYPE_CHECKING:
    from amici.sim.sundials.petab import PetabSimulationResult, PetabSimulator

__all__ = [
    "RobustConsistency",
    "run_simulation_to_cached_functions",
    "simulate_petab_to_cached_functions",
    "simulate_petab_v2_to_cached_functions",
]

LOG_E_10 = np.log(10)


class RobustConsistency(Consistency):
    """`Consistency`, plus rejection of step sizes that are self-consistent
    but inconsistent with the majority of other step sizes.

    `Consistency` checks whether the requested methods (e.g.
    forward/backward/central) agree with each other at each step size
    ("self-consistent"), then blends every self-consistent size's mean into
    the final value. Self-consistency alone is not a strong guarantee on its
    own: a step size can be small enough that all methods sample points
    within the target function's floating-point noise floor and become
    correlated (affected by the same rounding/cancellation error) --
    self-consistent, yet biased away from the truth. Symmetrically, a step
    size can also be large enough that all methods are biased the same way
    by higher-order/truncation effects.

    To guard against this, self-consistent step sizes are additionally
    required to agree with the majority of other self-consistent step sizes,
    via iterative outlier rejection (order-independent; step size magnitude
    is not used as a proxy for trustworthiness): repeatedly compute the
    median and a robust (MAD-based) spread of the current candidates, and
    drop the single worst-deviating one if it exceeds ``trend_n_sigma``
    scaled MADs from the median, until nothing looks anomalous. This only
    activates once there are at least ``min_trend_samples`` self-consistent
    step sizes; below that, there isn't enough data to estimate a spread, and
    all self-consistent step sizes are used, as in `Consistency`. A
    `UserWarning` is emitted whenever one or more step sizes are rejected
    this way.

    This addresses a long-standing intermittent CI failure in AMICI's PEtab
    benchmark gradient test
    (``test_benchmark_gradient[Weber_BMC2015-*-unscaled]``, see
    https://github.com/AMICI-dev/AMICI/issues/3078): that test uses
    `Consistency` to finite-difference-check an analytically computed
    gradient for a model parameter (``a32``) several orders of magnitude
    smaller than the model's other free parameters, and a small step size
    could become spuriously self-consistent while biased away from the true
    derivative.

    Note that this is a majority-vote style method: like any check based
    purely on the agreement of the values themselves (no independent ground
    truth), it has a breakdown point of roughly 50% (a property of the
    underlying median/MAD statistics) -- if close to half (or more) of the
    self-consistent step sizes are corrupted, this check cannot reliably
    tell which subset is trustworthy. This is a fundamental limitation of
    any purely data-driven consistency check, not something this
    implementation can detect or work around; sufficient step sizes with a
    real chance of being individually trustworthy should be provided.

    This was originally proposed upstream, in fiddy, as
    https://github.com/ICB-DCM/fiddy/pull/77, but was not merged; it lives
    here instead.
    """

    id = "robust_consistency"

    def __init__(
        self,
        *args,
        trend_n_sigma: float = 5.0,
        min_trend_samples: int = 3,
        **kwargs,
    ):
        """Construct.

        :param trend_n_sigma:
            The number of scaled median-absolute-deviations a
            self-consistent step size's estimate may deviate from the
            median of the other trusted step sizes' estimates, before it
            is rejected as an outlier.
        :param min_trend_samples:
            The minimum number of self-consistent step sizes required
            before the cross-step-size outlier rejection is attempted.
            Below this, all self-consistent step sizes are trusted, same
            as in `Consistency`.
        :param args:
            Positional arguments passed to `Consistency.__init__`.
        :param kwargs:
            Keyword arguments passed to `Consistency.__init__`
            (e.g. ``rtol``, ``atol``, ``equal_nan``).
        """
        super().__init__(*args, **kwargs)
        self.trend_n_sigma = trend_n_sigma
        self.min_trend_samples = min_trend_samples

    def _self_consistent_means(
        self, directional_derivative: DirectionalDerivative
    ) -> list[Type.DIRECTIONAL_DERIVATIVE]:
        """Group results by step size, and return the per-size mean for
        every step size whose requested methods agree with each other
        ("self-consistent") within ``rtol/2``, ``atol/2``."""
        computer_results = directional_derivative.get_computer_results()
        analysis_results = directional_derivative.get_analysis_results()
        results_by_size = {}
        for result in [*computer_results, *analysis_results]:
            size = result.metadata.get("size_absolute", None)
            if size is None:
                continue
            if size not in results_by_size:
                results_by_size[size] = {}
            if result.method_id in results_by_size[size]:
                raise ValueError(
                    f"Duplicate, and possibly conflicting, results for method "
                    f'"{result.method_id}" and size "{size}".',
                )
            results_by_size[size][result.method_id] = result.value

        self_consistent_means = []
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", "Mean of empty slice", RuntimeWarning
            )
            for results in results_by_size.values():
                values = list(results.values())
                mean = np.nanmean(values, axis=0)
                is_self_consistent = np.isclose(
                    values,
                    mean,
                    rtol=self.rtol / 2,
                    atol=self.atol / 2,
                    equal_nan=self.equal_nan,
                ).all()
                if is_self_consistent:
                    self_consistent_means.append(mean)
        return self_consistent_means

    def method(
        self, directional_derivative: DirectionalDerivative
    ) -> tuple[bool, float]:
        self_consistent_means = self._self_consistent_means(
            directional_derivative
        )

        if not self_consistent_means:
            return False, np.nan

        trusted_means = self._reject_outliers(self_consistent_means)

        if not trusted_means:
            return False, np.nan

        n_rejected = len(self_consistent_means) - len(trusted_means)
        if n_rejected:
            warnings.warn(
                f"{n_rejected} step size(s) were self-consistent (the "
                "requested methods agreed with each other) but were "
                "rejected as inconsistent with the majority of other step "
                "sizes; see `RobustConsistency`'s docstring.",
                stacklevel=2,
            )

        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", "Mean of empty slice", RuntimeWarning
            )
            value = np.nanmean(trusted_means, axis=0)

        success = (
            np.isclose(
                trusted_means,
                value,
                rtol=self.rtol,
                atol=self.atol,
                equal_nan=self.equal_nan,
            ).all()
            and not np.isnan(trusted_means).all()
        )
        return success, value

    def _reject_outliers(
        self, means: list[Type.DIRECTIONAL_DERIVATIVE]
    ) -> list[Type.DIRECTIONAL_DERIVATIVE]:
        """Iteratively reject step sizes whose estimate is an outlier.

        See the class docstring for the rationale. Order-independent: does
        not assume larger (or smaller) step sizes are inherently more
        trustworthy.

        :param means:
            The per-step-size mean estimates that passed the
            within-step-size self-consistency check.
        :return:
            The subset of `means` that are also mutually consistent with
            each other.
        """
        trusted = list(means)
        if len(trusted) < self.min_trend_samples:
            return trusted

        floor = max(self.atol / 2, np.finfo(float).tiny)
        while len(trusted) >= self.min_trend_samples:
            stacked = np.asarray(trusted, dtype=float)
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", "All-NaN", RuntimeWarning)
                center = np.nanmedian(stacked, axis=0)
                mad = np.nanmedian(np.abs(stacked - center), axis=0)
                scale = np.maximum(mad * 1.4826, floor)
                # One badness score per candidate, reduced across all
                # output dimensions (a candidate is an outlier if it
                # deviates too much in *any* output element).
                badness = np.nanmax(
                    (np.abs(stacked - center) / scale).reshape(
                        len(trusted), -1
                    ),
                    axis=1,
                )
                worst = int(np.nanargmax(badness))
            if badness[worst] > self.trend_n_sigma:
                trusted.pop(worst)
            else:
                break
        return trusted


def _transform_gradient_lin_to_lin(gradient_value, _):
    return gradient_value


def _transform_gradient_lin_to_log(gradient_value, parameter_value):
    return gradient_value * parameter_value


def _transform_gradient_lin_to_log10(gradient_value, parameter_value):
    return gradient_value * (parameter_value * LOG_E_10)


transforms = {
    LIN: _transform_gradient_lin_to_lin,
    LOG: _transform_gradient_lin_to_log,
    LOG10: _transform_gradient_lin_to_log10,
}


all_rdata_derivatives = {
    "x": "sx",
    "x0": "sx0",
    "x_ss": "sx_ss",
    "y": "sy",
    "sigmay": "ssigmay",
    "z": "sz",
    "rz": "srz",
    "sigmaz": "ssigmaz",
    "llh": "sllh",
    "sllh": "s2llh",
    "res": "sres",
}

# The dimension of the AMICI ReturnData that contains parameters.
# Should be shifted to the last dimension to be compatible with fiddy.
derivative_parameter_dimension = {
    "sx": 1,
    "sx0": 0,
    "sx_ss": 0,
    "sy": 1,
    "ssigmay": 1,
    # 'sz'      : ???,
    "srz": 2,
    # 'ssigmaz' : ???,
    "sllh": 0,
    "s2llh": 1,
    "sres": 1,
}


def _rdata_array_transpose(array: np.ndarray, variable: str) -> tuple[int]:
    if array.size == 0:
        return array
    original_parameter_dimension = derivative_parameter_dimension[variable]
    return np.moveaxis(array, original_parameter_dimension, -1)


default_derivatives = {
    k: v
    for k, v in all_rdata_derivatives.items()
    if v not in ["sz", "srz", "ssigmaz", "s2llh"]
}


def run_simulation_to_cached_functions(
    amici_model: AmiciModel,
    *,
    cache: bool = True,
    free_parameter_ids: list[str] = None,
    amici_solver: AmiciSolver = None,
    amici_edata: AmiciExpData = None,
    derivative_variables: list[str] = None,
):
    """Convert `run_simulation` to fiddy functions.

    :param amici_model:
        The AMICI model to simulate.
    :param amici_solver:
        The AMICI solver to use. If `None`, a new solver will be created from
        the model.
    :param amici_edata:
        The AMICI ExpData to use. If `None`, no data will be used.
    :param derivative_variables:
        The variables that derivatives will be computed or approximated for.
        See the keys of `all_rdata_derivatives` for options.
    :param free_parameter_ids:
        The IDs that correspond to the values in the free parameter vector that is
        simulated.
    :param cache:
        Whether to cache the function calls.
    :returns: function, derivatives and structure
    """
    if amici_solver is None:
        amici_solver = amici_model.create_solver()
    if free_parameter_ids is None:
        free_parameter_ids = amici_model.get_free_parameter_ids()
    if amici_edata is not None and amici_edata.free_parameters is not None:
        raise NotImplementedError(
            "Customization of parameter values inside AMICI ExpData."
        )
    chosen_derivatives = default_derivatives
    if derivative_variables is not None:
        chosen_derivatives = {
            k: all_rdata_derivatives[k] for k in derivative_variables
        }

    def run_amici_simulation(
        point: Type.POINT, order: SensitivityOrder
    ) -> ReturnData:
        problem_parameters = dict(zip(free_parameter_ids, point, strict=True))
        amici_model.set_free_parameter_by_id(problem_parameters)
        amici_solver.set_sensitivity_order(order)
        rdata = run_simulation(
            model=amici_model, solver=amici_solver, edata=amici_edata
        )
        return rdata

    def function(point: Type.POINT):
        rdata = run_amici_simulation(point=point, order=SensitivityOrder.none)
        outputs = {
            variable: fiddy_array(getattr(rdata, variable))
            for variable in chosen_derivatives
        }
        rdata_flat = np.concatenate(
            [output.flat for output in outputs.values()]
        )
        return rdata_flat

    def derivative(point: Type.POINT, return_dict: bool = False):
        rdata = run_amici_simulation(point=point, order=SensitivityOrder.first)
        outputs = {
            variable: _rdata_array_transpose(
                array=fiddy_array(getattr(rdata, derivative_variable)),
                variable=derivative_variable,
            )
            for variable, derivative_variable in chosen_derivatives.items()
        }
        rdata_flat = np.concatenate(
            [
                output_array.reshape(-1, output_array.shape[-1])
                for output_array in outputs.values()
            ],
            axis=0,
        )
        if return_dict:
            return outputs
        return rdata_flat

    if cache:
        function = CachedFunction(function)
        derivative = CachedFunction(derivative)

    # Get structure
    dummy_point = fiddy_array(
        [
            amici_model.get_free_parameter_by_id(par_id)
            for par_id in free_parameter_ids
        ]
    )
    dummy_rdata = run_amici_simulation(
        point=dummy_point, order=SensitivityOrder.first
    )

    structures = {
        "function": {variable: None for variable in chosen_derivatives},
        "derivative": {variable: None for variable in chosen_derivatives},
    }
    function_position = 0
    derivative_position = 0
    for variable, derivative_variable in chosen_derivatives.items():
        function_array = fiddy_array(getattr(dummy_rdata, variable))
        derivative_array = fiddy_array(
            getattr(dummy_rdata, derivative_variable)
        )
        structures["function"][variable] = (
            function_position,
            function_position + function_array.size,
            function_array.shape,
        )
        structures["derivative"][variable] = (
            derivative_position,
            derivative_position + derivative_array.size,
            derivative_array.shape,
        )
        function_position += function_array.size
        derivative_position += derivative_array.size

    return function, derivative, structures


# (start, stop, shape)
TYPE_STRUCTURE = tuple[int, int, tuple[int, ...]]


def simulate_petab_to_cached_functions(
    petab_problem: petab.Problem,
    *,
    amici_model: Model,
    free_parameter_ids: list[str] = None,
    cache: bool = True,
    precreate_edatas: bool = True,
    precreate_parameter_mapping: bool = True,
    simulate_petab: Callable[[Any], str] = None,
    **kwargs,
) -> tuple[Type.FUNCTION, Type.FUNCTION]:
    """
    Convert :func:`amici.sim.sundials.petab.v1.simulate_petab`
    (PEtab v1 simulations) to fiddy functions.

    Note that all gradients are provided on linear scale. The correction from
    `'log10'` scale is automatically done.

    :param amici_model:
        The AMICI model to simulate.
    :param simulate_petab:
        A method to simulate PEtab problems with AMICI, e.g.
        `amici.petab_objective.simulate_petab`.
    :param free_parameter_ids:
        The IDs of the parameters, in the order that parameter values will
        be supplied. Defaults to `petab_problem.parameter_df.index`.
    :param petab_problem:
        The PEtab problem.
    :param cache:
        Whether to cache the function call.
    :param precreate_edatas:
        Whether to create the AMICI measurements object in advance, to save
        time.
    :param precreate_parameter_mapping:
        Whether to create the AMICI parameter mapping object in advance, to
        save time.
    :param kwargs:
        Passed to `simulate_petab`.
    :returns:
        A tuple of:

        * 1: A method to compute the function at a point.
        * 2: A method to compute the gradient at a point.
    """
    if free_parameter_ids is None:
        free_parameter_ids = list(petab_problem.parameter_df.index)

    if simulate_petab is None:
        from amici.sim.sundials.petab.v1._simulations import simulate_petab

    edatas = None
    if precreate_edatas:
        edatas = create_edatas(
            amici_model=amici_model,
            petab_problem=petab_problem,
            simulation_conditions=petab_problem.get_simulation_conditions_from_measurement_df(),
        )

    parameter_mapping = None
    if precreate_parameter_mapping:
        parameter_mapping = create_parameter_mapping(
            petab_problem=petab_problem,
            simulation_conditions=petab_problem.get_simulation_conditions_from_measurement_df(),
            scaled_parameters=kwargs.get(
                "scaled_parameters",
                (
                    signature(simulate_petab)
                    .parameters["scaled_parameters"]
                    .default
                ),
            ),
            amici_model=amici_model,
        )

    precreated_kwargs = {
        "edatas": edatas,
        "parameter_mapping": parameter_mapping,
        "petab_problem": petab_problem,
    }
    precreated_kwargs = {
        k: v for k, v in precreated_kwargs.items() if v is not None
    }

    amici_solver = kwargs.pop("solver", amici_model.create_solver())

    simulate_petab_partial = partial(
        simulate_petab,
        amici_model=amici_model,
        **precreated_kwargs,
        **kwargs,
    )

    def simulate_petab_full(point: Type.POINT, order: SensitivityOrder):
        problem_parameters = dict(zip(free_parameter_ids, point, strict=True))
        amici_solver.set_sensitivity_order(order)
        result = simulate_petab_partial(
            problem_parameters=problem_parameters,
            solver=amici_solver,
        )
        return result

    def function(point: Type.POINT):
        output = simulate_petab_full(point, order=SensitivityOrder.none)
        result = output[LLH]
        return np.array(result)

    def derivative(point: Type.POINT) -> Type.POINT:
        result = simulate_petab_full(point, order=SensitivityOrder.first)

        if result[SLLH] is None:
            raise RuntimeError("Simulation failed.")

        sllh = np.array(
            [result[SLLH][parameter_id] for parameter_id in free_parameter_ids]
        )
        return sllh

    if cache:
        function = CachedFunction(function)
        derivative = CachedFunction(derivative)

    return function, derivative


def simulate_petab_v2_to_cached_functions(
    petab_simulator: PetabSimulator,
    *,
    free_parameter_ids: list[str] = None,
    cache: bool = True,
) -> tuple[Type.FUNCTION, Type.FUNCTION]:
    r"""Create fiddy functions for PetabSimulator.

    :param petab_simulator:
        The PEtab simulator to use.
    :param free_parameter_ids:
        The IDs of the parameters, in the order that parameter values will
        be supplied. Defaults to the estimated parameters of the PEtab problem.
    :param cache:
        Whether to cache the function call.
    :returns:
        tuple of:

        * 1: A method to compute the function at a point.
        * 2: A method to compute the gradient at a point.
    """
    if free_parameter_ids is None:
        free_parameter_ids = list(petab_simulator._petab_problem.x_free_ids)

    def simulate(
        point: Type.POINT, order: SensitivityOrder
    ) -> PetabSimulationResult:
        problem_parameters = dict(zip(free_parameter_ids, point, strict=True))
        petab_simulator.solver.set_sensitivity_order(order)

        result = petab_simulator.simulate(
            problem_parameters=problem_parameters,
        )
        return result

    def function(point: Type.POINT) -> np.ndarray:
        output = simulate(point, order=SensitivityOrder.none)
        result = output.llh
        return np.array(result)

    def derivative(point: Type.POINT) -> Type.POINT:
        result = simulate(point, order=SensitivityOrder.first)

        if result.sllh is None:
            raise RuntimeError("Simulation failed.")

        sllh = np.array(
            [result.sllh[parameter_id] for parameter_id in free_parameter_ids]
        )
        return sllh

    if cache:
        function = CachedFunction(function)
        derivative = CachedFunction(derivative)

    return function, derivative
