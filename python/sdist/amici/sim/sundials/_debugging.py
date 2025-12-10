"""Functions for debugging AMICI simulation failures."""

from __future__ import annotations

import numpy as np

import amici.sim.sundials


def get_model_for_preeq(
    model: amici.sim.sundials.Model, edata: amici.sim.sundials.ExpData
):
    """Get a model set-up to simulate the preequilibration condition as
    specified in `edata`.

    Useful for analyzing simulation errors during preequilibration.
    Simulating the returned model will reproduce the behavior of
    simulation-based preequilibration.

    Note that for models with events, the simulation results may differ.
    During preequilibration, event-handling is disabled. However, when
    simulating the returned model, event-handling will be enabled.
    For events triggered at fixed timepoints, this can be avoided by setting
    :meth:`t0 <amici.Model.set_t0>` to a timepoints after the last trigger
    timepoint.

    :param model:
        The model for which *edata* was generated.
    :param edata:
        The experimental data object with a preequilibration condition
        specified.
    :return:
        A copy of *model* with the same parameters, initial states, ...
        as *amici_model* for the preequilibration condition.
        Output timepoints are set to ``[inf]`` and will have to be adjusted.
    """
    model = model.clone()
    model.set_timepoints([np.inf])
    model.set_fixed_parameters(edata.fixed_parameters_pre_equilibration)
    if edata.pscale:
        model.set_parameter_scale(edata.pscale)
    if edata.free_parameters:
        model.set_free_parameters(edata.free_parameters)
    if edata.plist:
        model.set_parameter_list(edata.plist)
    model.set_initial_state(edata.x0)
    # has to be set *after* parameter list/scale!
    model.set_initial_state_sensitivities(edata.sx0)
    model.set_t0(edata.t_start)

    return model
