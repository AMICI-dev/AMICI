#include "amici/simulation_parameters.h"
#include "amici/misc.h"

#include <numeric>

namespace amici {

bool operator==(SimulationParameters const& a, SimulationParameters const& b) {
    return is_equal(a.fixed_parameters, b.fixed_parameters)
           && is_equal(
               a.fixed_parameters_pre_equilibration,
               b.fixed_parameters_pre_equilibration
           )
           && is_equal(
               a.fixed_parameters_presimulation,
               b.fixed_parameters_presimulation
           )
           && is_equal(a.free_parameters, b.free_parameters) && (a.plist == b.plist)
           && (a.pscale == b.pscale)
           && (a.reinitialize_fixed_parameter_initial_states
               == b.reinitialize_fixed_parameter_initial_states)
           && is_equal(a.sx0, b.sx0) && (a.t_presim == b.t_presim)
           && (a.t_start == b.t_start) && (a.timepoints == b.timepoints)
           && ((a.t_start_preeq == b.t_start_preeq)
               || (std::isnan(a.t_start_preeq) && std::isnan(b.t_start_preeq)));
}

void SimulationParameters::
    reinitialize_all_fixed_parameter_dependent_initial_states_for_presimulation(
        int const nx_rdata
    ) {
    reinitialization_state_idxs_presim.resize(nx_rdata);
    std::iota(
        reinitialization_state_idxs_presim.begin(),
        reinitialization_state_idxs_presim.end(), 0
    );
}

void SimulationParameters::
    reinitialize_all_fixed_parameter_dependent_initial_states_for_simulation(
        int const nx_rdata
    ) {
    reinitialization_state_idxs_sim.resize(nx_rdata);
    std::iota(
        reinitialization_state_idxs_sim.begin(),
        reinitialization_state_idxs_sim.end(), 0
    );
}

void SimulationParameters::
    reinitialize_all_fixed_parameter_dependent_initial_states(
        int const nx_rdata
    ) {
    reinitialize_all_fixed_parameter_dependent_initial_states_for_presimulation(
        nx_rdata
    );
    reinitialize_all_fixed_parameter_dependent_initial_states_for_simulation(
        nx_rdata
    );
}

} // namespace amici
