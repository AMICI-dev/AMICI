"""
Plotting
--------
Plotting related functions
"""
from . import ReturnDataView, Model

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from typing import Optional, Iterable

def plotStateTrajectories(
        rdata: ReturnDataView,
        state_indices: Optional[Iterable[int]] = None,
        ax: Optional[Axes] = None,
        model: Model = None
) -> None:
    """
    Plot state trajectories

    :param rdata:
        AMICI simulation results as returned by
        :func:`amici.amici.runAmiciSimulation`

    :param state_indices:
        Indices of states for which trajectories are to be plotted

    :param ax:
        matplotlib Axes instance to plot into
    
    :param model:
        amici model instance
    """
    if not ax:
        fig, ax = plt.subplots()
    if not state_indices:
        state_indices = range(rdata['x'].shape[1])
    for ix in state_indices:
        if model is None:
            label = f'x_{ix}'
        elif model.getStateNames()[ix] != '':
            label = model.getStateNames()[ix]
        else:
            label = model.getStateIds()[ix]
        ax.plot(rdata['t'], rdata['x'][:, ix], label = label)
        ax.set_xlabel('$t$ (s)')
        ax.set_ylabel('$x_i(t)$ (mmol/ml)')
        ax.legend()
        ax.set_title('State trajectories')
    
    
def plotObservableTrajectories(
        rdata: ReturnDataView,
        observable_indices: Optional[Iterable[int]] = None,
        ax: Optional[Axes] = None,
        model: Model = None
) -> None:
    """
    Plot observable trajectories

    :param rdata:
        AMICI simulation results as returned by
        :func:`amici.amici.runAmiciSimulation`

    :param observable_indices:
        Indices of observables for which trajectories are to be plotted

    :param ax:
        matplotlib Axes instance to plot into

    :param model:
        amici model instance
    """
    if not ax:
        fig, ax = plt.subplots()
    if not observable_indices:
        observable_indices = range(rdata['y'].shape[1])
    for iy in observable_indices:
        if model is None:
            label = f'y_{iy}'
        elif model.getObservableNames()[iy] != '':
            label = model.getObservableNames()[iy]
        else:
            label = model.getObservableIds()[iy]
        ax.plot(rdata['t'], rdata['y'][:, iy], label = label)
        ax.set_xlabel('$t$ (s)')
        ax.set_ylabel('$y_i(t)$ (AU)')
        ax.legend()
        ax.set_title('Observables')
