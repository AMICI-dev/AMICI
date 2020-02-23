"""
Plotting
--------
Plotting related functions
"""
from . import ReturnDataView

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from typing import Optional, Iterable


def plotStateTrajectories(
        rdata: ReturnDataView,
        state_indices: Optional[Iterable[int]] = None,
        ax: Optional[Axes] = None
) -> None:
    """
    Plot state trajectories

    :param rdata:
        AMICI simulation results as returned by amici.getSimulationResults()

    :param state_indices:
        Indices of states for which trajectories are to be plotted

    :param ax:
        matplotlib.axes.Axes instance to plot into

    """
    if not ax:
        fig, ax = plt.subplots()
    if not state_indices:
        state_indices = range(rdata['x'].shape[1])
    for ix in state_indices:
        ax.plot(rdata['t'], rdata['x'][:, ix], label='$x_{%d}$' % ix)
        ax.set_xlabel('$t$ (s)')
        ax.set_ylabel('$x_i(t)$ (mmol/ml)')
        ax.legend()
        ax.set_title('State trajectories')
    
    
def plotObservableTrajectories(
        rdata: ReturnDataView,
        observable_indices: Optional[Iterable[int]] = None,
        ax: Optional[Axes] = None
) -> None:
    """
    Plot observable trajectories

    :param rdata:
        AMICI simulation results as returned by amici.getSimulationResults()

    :param observable_indices:
        Indices of observables for which trajectories are to be plotted

    :param ax:
        matplotlib.axes.Axes instance to plot into

    """
    if not ax:
        fig, ax = plt.subplots()
    if not observable_indices:
        observable_indices = range(rdata['y'].shape[1])
    for iy in observable_indices:
        ax.plot(rdata['t'], rdata['y'][:, iy], label='$y_{%d}$' % iy)
        ax.set_xlabel('$t$ (s)')
        ax.set_ylabel('$y_i(t)$ (AU)')
        ax.legend()
        ax.set_title('Observables')
