"""@package amici.plotting Plotting related functions"""

import matplotlib.pyplot as plt 

def plotStateTrajectories(rdata, state_indices=None, ax = None, model = None):
    """Plot state trajectories
    
    Arguments:
    rdata: AMICI simulation results as returned by amici.getSimulationResults()
    state_indices: Indices of states for which trajectories are to be plotted
    ax: matplotlib.axes.Axes instance to plot into
    model: Model instance
    
    Returns:

    Raises:

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
    
    
def plotObservableTrajectories(rdata, observable_indices=None, ax = None, model = None):
    """Plot observable trajectories
    
    Arguments:
    rdata: AMICI simulation results as returned by amici.getSimulationResults()
    observable_indices: Indices of observables for which trajectories are to be plotted
    ax: matplotlib.axes.Axes instance to plot into
    model: Model instance

    Returns:

    Raises:

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
