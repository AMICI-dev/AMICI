#include "amici/simulation_parameters.h"

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


} // namespace amici
