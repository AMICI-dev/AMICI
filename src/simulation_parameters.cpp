#include "amici/simulation_parameters.h"

#include <numeric>

namespace amici {

bool operator==(const SimulationParameters &a, const SimulationParameters &b) {
    return (a.fixedParameters == b.fixedParameters) &&
            (a.fixedParametersPreequilibration == b.fixedParametersPreequilibration) &&
            (a.fixedParametersPresimulation == b.fixedParametersPresimulation) &&
            (a.parameters == b.parameters) &&
            (a.plist == b.plist) &&
            (a.pscale == b.pscale) &&
            (a.reinitializeFixedParameterInitialStates == b.reinitializeFixedParameterInitialStates) &&
            (a.sx0 == b.sx0) &&
            (a.t_presim == b.t_presim) &&
            (a.tstart_ == b.tstart_) &&
            (a.ts_ == b.ts_);
}

void SimulationParameters::reinitializeAllFixedParameterDependentInitialStatesForPresimulation(int nx_rdata)
{
    reinitialization_state_idxs_presim.resize(nx_rdata);
    std::iota(reinitialization_state_idxs_presim.begin(),
              reinitialization_state_idxs_presim.end(), 0);

}

void SimulationParameters::reinitializeAllFixedParameterDependentInitialStatesForSimulation(int nx_rdata)
{
    reinitialization_state_idxs_sim.resize(nx_rdata);
    std::iota(reinitialization_state_idxs_sim.begin(),
              reinitialization_state_idxs_sim.end(), 0);
}

void SimulationParameters::reinitializeAllFixedParameterDependentInitialStates(int nx_rdata)
{
    reinitializeAllFixedParameterDependentInitialStatesForPresimulation(nx_rdata);
    reinitializeAllFixedParameterDependentInitialStatesForSimulation(nx_rdata);
}


} // namespace amici
