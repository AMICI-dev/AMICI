#ifndef AMICI_SIMULATION_PARAMETERS_H
#define AMICI_SIMULATION_PARAMETERS_H

#include "amici/defines.h"

#include <vector>

namespace amici {

/**
 * @brief Container for various simulation parameters.
 */
class SimulationParameters {
  public:
    SimulationParameters() = default;

    /**
     * @brief Constructor
     *
     * @param timepoints Timepoints for which simulation results are requested
     */
    explicit SimulationParameters(std::vector<realtype> timepoints)
        : timepoints(std::move(timepoints)) {}

    /**
     * @brief Constructor
     *
     * @param fixed_parameters Model parameters excluded from sensitivity
     * analysis
     * @param free_parameters Model parameters included in sensitivity analysis
     */
    SimulationParameters(
        std::vector<realtype> fixed_parameters,
        std::vector<realtype> free_parameters
    )
        : fixed_parameters(std::move(fixed_parameters))
        , free_parameters(std::move(free_parameters))
        , pscale(
              std::vector(this->free_parameters.size(), ParameterScaling::none)
          ) {}

#ifndef SWIGPYTHON
    /*
     * include/amici/simulation_parameters.h:71: Warning 509: Overloaded method
     * amici::SimulationParameters::SimulationParameters(std::vector<
     * amici::realtype,std::allocator< amici::realtype > >,std::vector<
     * amici::realtype,std::allocator< amici::realtype > >,std::vector<
     * amici::realtype,std::allocator< amici::realtype > >) effectively ignored,
     * include/amici/simulation_parameters.h:54: Warning 509: as it is shadowed
     * by amici::SimulationParameters::SimulationParameters(std::vector<
     * amici::realtype,std::allocator< amici::realtype > >,std::vector<
     * amici::realtype,std::allocator< amici::realtype > >,std::vector<
     * int,std::allocator< int > >).
     */
    /**
     * @brief Constructor
     *
     * @param fixed_parameters Model parameters excluded from sensitivity
     * analysis
     * @param free_parameters Model parameters included in sensitivity analysis
     * @param plist Model parameter indices w.r.t. which sensitivities are to be
     * computed
     */
    SimulationParameters(
        std::vector<realtype> fixed_parameters,
        std::vector<realtype> free_parameters, std::vector<int> plist
    )
        : fixed_parameters(std::move(fixed_parameters))
        , free_parameters(std::move(free_parameters))
        , pscale(
              std::vector(this->free_parameters.size(), ParameterScaling::none)
          )
        , plist(std::move(plist)) {}

    /**
     * @brief Constructor
     *
     * @param timepoints Timepoints for which simulation results are requested
     * @param fixed_parameters Model parameters excluded from sensitivity
     * analysis
     * @param free_parameters Model parameters included in sensitivity analysis
     */
    SimulationParameters(
        std::vector<realtype> timepoints,
        std::vector<realtype> fixed_parameters,
        std::vector<realtype> free_parameters
    )
        : fixed_parameters(std::move(fixed_parameters))
        , free_parameters(std::move(free_parameters))
        , pscale(
              std::vector(this->free_parameters.size(), ParameterScaling::none)
          )
        , timepoints(std::move(timepoints)) {}
#endif

    /**
     * @brief Set reinitialization of all states based on model constants for
     * presimulation (only meaningful if preequilibration is performed).
     *
     * Convenience function to populate
     * `reinitialization_state_idxs_presim` and
     * `reinitialization_state_idxs_sim`
     *
     * @param nx_rdata Number of states (Model::nx_rdata)
     */
    void
    reinitialize_all_fixed_parameter_dependent_initial_states_for_presimulation(
        int nx_rdata
    );

    /**
     * @brief Set reinitialization of all states based on model constants for
     * the 'main' simulation (only meaningful if presimulation or
     * preequilibration is performed).
     *
     * Convenience function to populate
     * `reinitialization_state_idxs_presim` and
     * `reinitialization_state_idxs_sim`
     *
     * @param nx_rdata Number of states (Model::nx_rdata)
     */
    void
    reinitialize_all_fixed_parameter_dependent_initial_states_for_simulation(
        int nx_rdata
    );

    /**
     * @brief Set reinitialization of all states based on model constants for
     * all simulation phases.
     *
     * Convenience function to populate
     * `reinitialization_state_idxs_presim` and
     * `reinitialization_state_idxs_sim`
     *
     * @param nx_rdata Number of states (Model::nx_rdata)
     */
    void
    reinitialize_all_fixed_parameter_dependent_initial_states(int nx_rdata);

    /**
     * @brief Model constants
     *
     * Vector of size Model::nk() or empty
     */
    std::vector<realtype> fixed_parameters;

    /**
     * @brief Model constants for pre-equilibration
     *
     * Vector of size Model::nk() or empty.
     */
    std::vector<realtype> fixed_parameters_pre_equilibration;

    /**
     * @brief Model constants for pre-simulation
     *
     * Vector of size Model::nk() or empty.
     */
    std::vector<realtype> fixed_parameters_presimulation;

    /**
     * @brief Model free_parameters
     *
     * Vector of size Model::np() or empty with parameter scaled according to
     * SimulationParameter::pscale.
     */
    std::vector<realtype> free_parameters;

    /**
     * @brief Initial state
     *
     * Vector of size Model::nx() or empty
     */
    std::vector<realtype> x0;

    /**
     * @brief Initial state sensitivities
     *
     * Dimensions:
     * Model::nx() * Model::nplist(), Model::nx() * ExpData::plist.size(), if
     * ExpData::plist is not empty, or empty
     */
    std::vector<realtype> sx0;

    /**
     * @brief Parameter scales
     *
     * Vector of parameter scale of size Model::np(), indicating how/if each
     * parameter is to be scaled.
     */
    std::vector<ParameterScaling> pscale;

    /**
     * @brief Parameter indices w.r.t. which to compute sensitivities
     */
    std::vector<int> plist;

    /**
     * @brief The initial time for pre-equilibration..
     *
     * NAN indicates that `tstart_` should be used.
     */
    realtype t_start_preeq{NAN};

    /**
     * @brief Starting time of the simulation.
     *
     * Output timepoints are absolute timepoints, independent of
     * \f$ t_{start} \f$.
     * For output timepoints \f$ t <  t_{start} \f$, the initial state will be
     * returned.
     */
    realtype t_start{0.0};

    /**
     * @brief Duration of pre-simulation.
     *
     * If this is > 0, presimulation will be performed from
     * (model->t0 - t_presim) to model->t0 using the fixed_parameters in
     * fixed_parameters_presimulation
     */
    realtype t_presim{0.0};

    /**
     * @brief Timepoints for which model state/outputs/... are requested
     *
     * Vector of timepoints.
     */
    std::vector<realtype> timepoints;

    /**
     * @brief Flag indicating whether reinitialization of states depending on
     * fixed parameters is activated
     */
    bool reinitialize_fixed_parameter_initial_states{false};

    /**
     * @brief Indices of states to be reinitialized based on provided
     * presimulation constants / fixed parameters.
     */
    std::vector<int> reinitialization_state_idxs_presim;

    /**
     * @brief Indices of states to be reinitialized based on provided
     * constants / fixed parameters.
     */
    std::vector<int> reinitialization_state_idxs_sim;
};

bool operator==(SimulationParameters const& a, SimulationParameters const& b);

} // namespace amici

#endif // AMICI_SIMULATION_PARAMETERS_H
