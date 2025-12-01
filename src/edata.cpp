#include "amici/edata.h"
#include "amici/defines.h"
#include "amici/model.h"
#include "amici/rdata.h"
#include "amici/symbolic_functions.h" // get_nan

#include <algorithm>
#include <random>
#include <utility>

namespace amici {

using std::isnan;

ExpData::ExpData(int const nytrue, int const nztrue, int const nmaxevent)
    : nytrue_(nytrue)
    , nztrue_(nztrue)
    , nmaxevent_(nmaxevent) {
    apply_dimensions();
}

ExpData::ExpData(
    int const nytrue, int const nztrue, int const nmaxevent,
    std::vector<realtype> ts
)
    : SimulationParameters(ts)
    , nytrue_(nytrue)
    , nztrue_(nztrue)
    , nmaxevent_(nmaxevent) {
    apply_dimensions();
}

ExpData::ExpData(
    int const nytrue, int const nztrue, int const nmaxevent,
    std::vector<realtype> ts, std::vector<realtype> fixed_parameters
)
    : SimulationParameters(ts)
    , nytrue_(nytrue)
    , nztrue_(nztrue)
    , nmaxevent_(nmaxevent) {
    this->fixed_parameters = std::move(fixed_parameters);
    apply_dimensions();
}

ExpData::ExpData(
    int const nytrue, int const nztrue, int const nmaxevent,
    std::vector<realtype> ts, std::vector<realtype> const& my,
    std::vector<realtype> const& sigma_y,
    std::vector<realtype> const& mz,
    std::vector<realtype> const& sigma_z
)
    : SimulationParameters(ts)
    , nytrue_(nytrue)
    , nztrue_(nztrue)
    , nmaxevent_(nmaxevent) {
    apply_dimensions();
    set_measurements(my);
    set_measurement_error(sigma_y);
    set_event_measurements(mz);
    set_event_measurement_error(sigma_z);
}

ExpData::ExpData(Model const& model)
    : ExpData(
          model.nytrue, model.nztrue, model.n_max_event(),
          model.get_timepoints(), model.get_fixed_parameters()
      ) {
    reinitialize_fixed_parameter_initial_states
        = model.get_reinitialize_fixed_parameter_initial_states()
          && model.get_reinitialization_state_idxs().empty();
    reinitialization_state_idxs_sim = model.get_reinitialization_state_idxs();
}

ExpData::ExpData(
    ReturnData const& rdata, realtype const sigma_y, realtype const sigma_z,
    int const seed
)
    : ExpData(
          rdata, std::vector<realtype>(rdata.nytrue * rdata.nt, sigma_y),
          std::vector<realtype>(rdata.nztrue * rdata.nmaxevent, sigma_z), seed
      ) {}

ExpData::ExpData(
    ReturnData const& rdata, std::vector<realtype> sigma_y,
    std::vector<realtype> sigma_z, int const seed
)
    : ExpData(rdata.nytrue, rdata.nztrue, rdata.nmaxevent, rdata.ts) {
    if (sigma_y.size() != static_cast<unsigned>(nytrue_)
        && sigma_y.size() != static_cast<unsigned>(nytrue_) * nt())
        throw AmiException(
            "Dimension of sigma_y must be %d or %d, was %d", nytrue_,
            nytrue_ * nt(), sigma_y.size()
        );

    if (sigma_z.size() != static_cast<unsigned>(nztrue_)
        && sigma_z.size() != static_cast<unsigned>(nztrue_) * nmaxevent_)
        throw AmiException(
            "Dimension of sigma_z must be %d or %d, was %d", nztrue_,
            nztrue_ * nmaxevent_, sigma_z.size()
        );

    std::mt19937 gen{seed < 0 ? std::random_device()() : seed};

    realtype sigma;

    check_sigma_positivity(sigma_y, "sigma_y");
    check_sigma_positivity(sigma_z, "sigma_z");

    for (int iy = 0; iy < nytrue_; ++iy) {
        for (int it = 0; it < nt(); ++it) {
            sigma = sigma_y.size() == static_cast<unsigned>(nytrue_)
                        ? sigma_y.at(iy)
                        : sigma_y.at(iy + nytrue_ * it);
            std::normal_distribution<> e{0, sigma};
            measurements_.at(iy + nytrue_ * it)
                = rdata.y.at(iy + rdata.ny * it) + e(gen);
            measurement_error_.at(iy + nytrue_ * it) = sigma;
        }
    }

    for (int iz = 0; iz < nztrue_; ++iz) {
        for (int ie = 0; ie < nmaxevent_; ++ie) {
            sigma = sigma_z.size() == static_cast<unsigned>(nztrue_)
                        ? sigma_z.at(iz)
                        : sigma_z.at(iz + nztrue_ * ie);
            std::normal_distribution<> e{0, sigma};
            event_measurements.at(iz + rdata.nztrue * ie)
                = rdata.z.at(iz + rdata.nz * ie) + e(gen);
            event_measurementserror_.at(iz + rdata.nztrue * ie) = sigma;
        }
    }

    id = rdata.id;
}

void ExpData::set_timepoints(std::vector<realtype> const& ts) {
    if (!std::ranges::is_sorted(ts))
        throw AmiException(
            "Encountered non-monotonic timepoints, please order timepoints "
            "such that they are monotonically increasing!"
        );
    timepoints = ts;
    apply_data_dimension();
}

std::vector<realtype> const& ExpData::get_timepoints() const {
    return timepoints;
}

int ExpData::nt() const { return gsl::narrow<int>(timepoints.size()); }

realtype ExpData::get_timepoint(int const it) const {
    return timepoints.at(it);
}

void ExpData::set_measurements(std::vector<realtype> const& measurements) {
    check_data_dimension(measurements, "measurements");

    if (measurements.size() == static_cast<unsigned>(nt()) * nytrue_)
        measurements_ = measurements;
    else if (measurements.empty())
        measurements_.clear();
}

void ExpData::set_measurements(
    std::vector<realtype> const& measurements, int iy
) {
    if (measurements.size() != static_cast<unsigned>(nt()))
        throw AmiException(
            "Input measurements did not match dimensions nt (%i), was %i", nt(),
            measurements.size()
        );

    for (int it = 0; it < nt(); ++it)
        measurements_.at(iy + it * nytrue_) = measurements.at(it);
}

bool ExpData::is_set_measurements(int const it, int const iy) const {
    return !measurements_.empty()
           && !isnan(measurements_.at(it * nytrue_ + iy));
}

std::vector<realtype> const& ExpData::get_measurements() const {
    return measurements_;
}

realtype const* ExpData::get_measurements_ptr(int const it) const {
    if (!measurements_.empty())
        return &measurements_.at(it * nytrue_);

    return nullptr;
}

void ExpData::set_measurement_error(
    std::vector<realtype> const& sigma_y
) {
    check_data_dimension(sigma_y, "sigma_y");
    check_sigma_positivity(sigma_y, "sigma_y");

    if (sigma_y.size() == static_cast<unsigned>(nt()) * nytrue_)
        measurement_error_ = sigma_y;
    else if (sigma_y.empty())
        measurement_error_.clear();
}

void ExpData::set_measurement_error(realtype const sigma) {
    check_sigma_positivity(sigma, "sigma");
    std::ranges::fill(measurement_error_, sigma);
}

void ExpData::set_measurement_error(
    std::vector<realtype> const& sigma_y, int const iy
) {
    if (sigma_y.size() != static_cast<unsigned>(nt()))
        throw AmiException(
            "Input sigma_y did not match dimensions nt (%i), was %i",
            nt(), sigma_y.size()
        );
    check_sigma_positivity(sigma_y, "sigma_y");

    for (int it = 0; it < nt(); ++it)
        measurement_error_.at(iy + it * nytrue_)
            = sigma_y.at(it);
}

void ExpData::set_measurement_error(realtype const sigma, int const iy) {
    check_sigma_positivity(sigma, "sigma");
    for (int it = 0; it < nt(); ++it)
        measurement_error_.at(iy + it * nytrue_) = sigma;
}

bool ExpData::is_set_measurement_error(int const it, int const iy) const {
    return !measurement_error_.empty()
           && !isnan(measurement_error_.at(it * nytrue_ + iy));
}

std::vector<realtype> const& ExpData::get_measurement_error() const {
    return measurement_error_;
}

realtype const* ExpData::get_measurement_error_ptr(int const it) const {
    if (!measurement_error_.empty())
        return &measurement_error_.at(it * nytrue_);

    return nullptr;
}

void ExpData::set_event_measurements(std::vector<realtype> const& event_measurements) {
    check_events_dimension(event_measurements, "event_measurements");

    if (event_measurements.size() == static_cast<unsigned>(nmaxevent_) * nztrue_)
        event_measurements_ = event_measurements;
    else if (event_measurements.empty())
        event_measurements_.clear();
}

void ExpData::set_event_measurements(
    std::vector<realtype> const& event_measurements, int const iz
) {
    if (event_measurements.size() != static_cast<unsigned>(nmaxevent_)) {
        throw AmiException(
            "Input event_measurements did not match dimensions nmaxevent (%i), was "
            "%i",
            nmaxevent_, event_measurements.size()
        );
    }

    for (int ie = 0; ie < nmaxevent_; ++ie)
        event_measurements_.at(iz + ie * nztrue_) = event_measurements.at(ie);
}

bool ExpData::is_set_event_measurements(int const ie, int const iz) const {
    return !event_measurements_.empty()
           && !isnan(event_measurements_.at(ie * nztrue_ + iz));
}

std::vector<realtype> const& ExpData::get_event_measurements() const {
    return event_measurements_;
}

realtype const* ExpData::get_event_measurements_ptr(int const ie) const {
    if (!event_measurements_.empty())
        return &event_measurements_.at(ie * nztrue_);

    return nullptr;
}

void ExpData::set_event_measurement_error(
    std::vector<realtype> const& sigma_z
) {
    check_events_dimension(sigma_z, "sigma_z");
    check_sigma_positivity(sigma_z, "sigma_z");

    if (sigma_z.size() == (unsigned)nmaxevent_ * nztrue_)
        event_measurement_error_ = sigma_z;
    else if (sigma_z.empty())
        event_measurement_error_.clear();
}

void ExpData::set_event_measurement_error(realtype const sigma) {
    check_sigma_positivity(sigma, "sigma");
    std::ranges::fill(event_measurement_error_, sigma);
}

void ExpData::set_event_measurement_error(
    std::vector<realtype> const& sigma_z, int const iz
) {
    if (sigma_z.size() != (unsigned)nmaxevent_)
        throw AmiException(
            "Input sigma_z did not match dimensions nmaxevent "
            "(%i), was %i",
            nmaxevent_, sigma_z.size()
        );
    check_sigma_positivity(sigma_z, "sigma_z");

    for (int ie = 0; ie < nmaxevent_; ++ie)
        event_measurement_error_.at(iz + ie * nztrue_)
            = sigma_z.at(ie);
}

void ExpData::set_event_measurement_error(realtype const sigma, int const iz) {
    check_sigma_positivity(sigma, "sigma");

    for (int ie = 0; ie < nmaxevent_; ++ie)
        event_measurement_error_.at(iz + ie * nztrue_) = sigma;
}

bool ExpData::is_set_event_measurement_error(int const ie, int const iz) const {
    if (!event_measurement_error_.empty()) // avoid out of bound memory access
        return !isnan(event_measurement_error_.at(ie * nztrue_ + iz));

    return false;
}

std::vector<realtype> const& ExpData::get_event_measurement_error() const {
    return event_measurement_error_;
}

realtype const* ExpData::get_event_measurement_error_ptr(int const ie) const {
    if (!event_measurement_error_.empty())
        return &event_measurement_error_.at(ie * nztrue_);

    return nullptr;
}

void ExpData::clear_observations() {
    std::ranges::fill(measurements_, get_nan());
    std::ranges::fill(measurement_error_, get_nan());
    std::ranges::fill(event_measurements_, get_nan());
    std::ranges::fill(event_measurement_error_, get_nan());
}

void ExpData::apply_dimensions() {
    apply_data_dimension();
    apply_event_dimension();
}

void ExpData::apply_data_dimension() {
    measurements_.resize(nt() * nytrue_, get_nan());
    measurement_error_.resize(nt() * nytrue_, get_nan());
}

void ExpData::apply_event_dimension() {
    event_measurements.resize(nmaxevent_ * nztrue_, get_nan());
    event_measurement_error_.resize(nmaxevent_ * nztrue_, get_nan());
}

void ExpData::check_data_dimension(
    std::vector<realtype> const& input, char const* fieldname
) const {
    if (input.size() != static_cast<unsigned>(nt()) * nytrue_ && !input.empty())
        throw AmiException(
            "Input %s did not match dimensions nt (%i) x nytrue (%i), was %i",
            fieldname, nt(), nytrue_, input.size()
        );
}

void ExpData::check_events_dimension(
    std::vector<realtype> const& input, char const* fieldname
) const {
    if (input.size() != static_cast<unsigned>(nmaxevent_) * nztrue_
        && !input.empty())
        throw AmiException(
            "Input %s did not match dimensions nt (%i) x nytrue (%i), was %i",
            fieldname, nmaxevent_, nztrue_, input.size()
        );
}

void check_sigma_positivity(
    std::vector<realtype> const& sigmaVector, char const* vectorName
) {
    for (auto&& sigma : sigmaVector)
        check_sigma_positivity(sigma, vectorName);
}

void check_sigma_positivity(realtype const sigma, char const* sigmaName) {
    if (sigma <= 0.0)
        throw AmiException(
            "Encountered sigma <= 0 in %s! value: %f", sigmaName, sigma
        );
}

int ExpData::nytrue() const { return nytrue_; }

int ExpData::nztrue() const { return nztrue_; }

int ExpData::nmaxevent() const { return nmaxevent_; }

ConditionContext::ConditionContext(
    Model* model, ExpData const* edata, FixedParameterContext fpc
)
    : model_(model)
    , original_parameters_(model->get_free_parameters())
    , original_fixed_parameters_(model->get_fixed_parameters())
    , original_tstart_(model->t0())
    , original_tstart_preeq_(model->t0_preeq())
    , original_timepoints_(model->get_timepoints())
    , original_parameter_list_(model->get_parameter_list())
    , original_scaling_(model->get_parameter_scale())
    , original_reinitialize_fixed_parameter_initial_states_(
          model->get_reinitialize_fixed_parameter_initial_states()
          && model->get_reinitialization_state_idxs().empty()
      )
    , original_reinitialization_state_idxs(
          model->get_reinitialization_state_idxs()
      ) {
    if (model->has_custom_initial_state())
        original_x0_ = model->get_initial_state();

    if (model->has_custom_initial_state_sensitivities())
        original_sx0_ = model->get_initial_state_sensitivities();

    apply_condition(edata, fpc);
}

ConditionContext::~ConditionContext() { restore(); }

void ConditionContext::apply_condition(
    ExpData const* edata, FixedParameterContext fpc
) {
    if (!edata)
        return;

    // this needs to go first, otherwise nplist will not have the right
    // dimension for all other fields that depend on Model::nplist
    if (!edata->plist.empty())
        model_->set_parameter_list(edata->plist);

    // this needs to go second as setParameterScale will reset sx0
    if (!edata->pscale.empty()) {
        if (edata->pscale.size() != (unsigned)model_->np())
            throw AmiException(
                "Number of parameters (%d) in model does not"
                " match ExpData (%zd).",
                model_->np(), edata->pscale.size()
            );
        model_->set_parameter_scale(edata->pscale);
    }

    // this needs to be set in the model before handling initial state
    // sensitivities, which may be unscaled using model parameter values
    if (!edata->free_parameters.empty()) {
        if (edata->free_parameters.size() != (unsigned)model_->np())
            throw AmiException(
                "Number of parameters (%d) in model does not"
                " match ExpData (%zd).",
                model_->np(), edata->free_parameters.size()
            );
        model_->set_free_parameters(edata->free_parameters);
    }

    if (!edata->x0.empty()) {
        if (edata->x0.size() != (unsigned)model_->nx_rdata)
            throw AmiException(
                "Number of initial conditions (%d) in model does"
                " not match ExpData (%zd).",
                model_->nx_rdata, edata->x0.size()
            );
        model_->set_initial_state(edata->x0);
    }

    if (!edata->sx0.empty()) {
        if (edata->sx0.size() != (unsigned)model_->nx_rdata * model_->nplist())
            throw AmiException(
                "Number of initial conditions sensitivities (%d)"
                " in model does not match ExpData (%zd).",
                model_->nx_rdata * model_->nplist(), edata->sx0.size()
            );
        model_->set_initial_state_sensitivities(edata->sx0);
    }

    model_->set_reinitialize_fixed_parameter_initial_states(
        edata->reinitialize_fixed_parameter_initial_states
    );

    switch (fpc) {
    case FixedParameterContext::simulation:
        if (!edata->fixed_parameters.empty()) {
            // fixed parameter in model are superseded by those provided in
            // edata
            if (edata->fixed_parameters.size() != (unsigned)model_->nk())
                throw AmiException(
                    "Number of fixed parameters (%d) in model does"
                    "not match ExpData (%zd).",
                    model_->nk(), edata->fixed_parameters.size()
                );
            model_->set_fixed_parameters(edata->fixed_parameters);
            if (!edata->reinitialize_fixed_parameter_initial_states)
                model_->set_reinitialization_state_idxs(
                    edata->reinitialization_state_idxs_sim
                );
        }
        break;
    case FixedParameterContext::preequilibration:
        if (!edata->fixed_parameters_pre_equilibration.empty()) {
            // fixed parameter in model are superseded by those provided in
            // edata
            if (edata->fixed_parameters_pre_equilibration.size()
                != (unsigned)model_->nk())
                throw AmiException(
                    "Number of fixed parameters (%d) in model does"
                    "not match ExpData (preequilibration) (%zd).",
                    model_->nk(),
                    edata->fixed_parameters_pre_equilibration.size()
                );
            model_->set_fixed_parameters(
                edata->fixed_parameters_pre_equilibration
            );
        }
        break;
    case FixedParameterContext::presimulation:
        if (!edata->fixed_parameters_presimulation.empty()) {
            // fixed parameter in model are superseded by those provided in
            // edata
            if (edata->fixed_parameters_presimulation.size()
                != (unsigned)model_->nk())
                throw AmiException(
                    "Number of fixed parameters (%d) in model does"
                    " not match ExpData (presimulation) (%zd).",
                    model_->nk(), edata->fixed_parameters_presimulation.size()
                );
            model_->set_fixed_parameters(edata->fixed_parameters_presimulation);
            if (!edata->reinitialize_fixed_parameter_initial_states)
                model_->set_reinitialization_state_idxs(
                    edata->reinitialization_state_idxs_presim
                );
        }
        break;
    }

    model_->set_t0(edata->t_start);
    if (edata->nt()) {
        // fixed parameter in model are superseded by those provided in edata
        model_->set_timepoints(edata->get_timepoints());
    }
}

void ConditionContext::restore() {
    // parameter list has to be set before initial state sensitivities
    model_->set_parameter_list(original_parameter_list_);
    // parameter scale has to be set before initial state sensitivities
    model_->set_parameter_scale(original_scaling_);

    if (!original_x0_.empty())
        model_->set_initial_state(original_x0_);

    if (!original_sx0_.empty())
        model_->set_unscaled_initial_state_sensitivities(original_sx0_);

    model_->set_free_parameters(original_parameters_);
    model_->set_fixed_parameters(original_fixed_parameters_);
    model_->set_t0(original_tstart_);
    model_->set_t0_preeq(original_tstart_preeq_);
    model_->set_timepoints(original_timepoints_);
    model_->set_reinitialize_fixed_parameter_initial_states(
        original_reinitialize_fixed_parameter_initial_states_
    );
    model_->set_reinitialization_state_idxs(
        original_reinitialization_state_idxs
    );
}

} // namespace amici
