#include "amici/simulation_parameters.h"
#include "amici/misc.h"

#include <numeric>

namespace amici {

bool operator==(SimulationParameters const& a, SimulationParameters const& b) {
    return is_equal(a.fixedParameters, b.fixedParameters)
           && is_equal(
               a.fixedParametersPreequilibration,
               b.fixedParametersPreequilibration
           )
           && is_equal(
               a.fixedParametersPresimulation, b.fixedParametersPresimulation
           )
           && is_equal(a.parameters, b.parameters) && (a.plist == b.plist)
           && (a.pscale == b.pscale)
           && (a.reinitializeFixedParameterInitialStates
               == b.reinitializeFixedParameterInitialStates)
           && is_equal(a.sx0, b.sx0) && (a.t_presim == b.t_presim)
           && (a.tstart_ == b.tstart_) && (a.ts_ == b.ts_);
}

void SimulationParameters::
    reinitializeAllFixedParameterDependentInitialStatesForPresimulation(
        int nx_rdata
    ) {
    reinitialization_state_idxs_presim.resize(nx_rdata);
    std::iota(
        reinitialization_state_idxs_presim.begin(),
        reinitialization_state_idxs_presim.end(), 0
    );
}

void SimulationParameters::
    reinitializeAllFixedParameterDependentInitialStatesForSimulation(
        int nx_rdata
    ) {
    reinitialization_state_idxs_sim.resize(nx_rdata);
    std::iota(
        reinitialization_state_idxs_sim.begin(),
        reinitialization_state_idxs_sim.end(), 0
    );
}

void SimulationParameters::reinitializeAllFixedParameterDependentInitialStates(
    int nx_rdata
) {
    reinitializeAllFixedParameterDependentInitialStatesForPresimulation(nx_rdata
    );
    reinitializeAllFixedParameterDependentInitialStatesForSimulation(nx_rdata);
}

} // namespace amici
