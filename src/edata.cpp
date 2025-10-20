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
    std::vector<realtype> ts, std::vector<realtype> fixedParameters
)
    : SimulationParameters(ts)
    , nytrue_(nytrue)
    , nztrue_(nztrue)
    , nmaxevent_(nmaxevent) {
    this->fixed_parameters = std::move(fixedParameters);
    apply_dimensions();
}

ExpData::ExpData(
    int const nytrue, int const nztrue, int const nmaxevent,
    std::vector<realtype> ts, std::vector<realtype> const& observedData,
    std::vector<realtype> const& observedDataStdDev,
    std::vector<realtype> const& observedEvents,
    std::vector<realtype> const& observedEventsStdDev
)
    : SimulationParameters(ts)
    , nytrue_(nytrue)
    , nztrue_(nztrue)
    , nmaxevent_(nmaxevent) {
    apply_dimensions();
    set_observed_data(observedData);
    set_observed_data_std_dev(observedDataStdDev);
    set_observed_events(observedEvents);
    set_observed_events_std_dev(observedEventsStdDev);
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

    checkSigmaPositivity(sigma_y, "sigma_y");
    checkSigmaPositivity(sigma_z, "sigma_z");

    for (int iy = 0; iy < nytrue_; ++iy) {
        for (int it = 0; it < nt(); ++it) {
            sigma = sigma_y.size() == static_cast<unsigned>(nytrue_)
                        ? sigma_y.at(iy)
                        : sigma_y.at(iy + nytrue_ * it);
            std::normal_distribution<> e{0, sigma};
            observed_data_.at(iy + nytrue_ * it)
                = rdata.y.at(iy + rdata.ny * it) + e(gen);
            observed_data_std_dev_.at(iy + nytrue_ * it) = sigma;
        }
    }

    for (int iz = 0; iz < nztrue_; ++iz) {
        for (int ie = 0; ie < nmaxevent_; ++ie) {
            sigma = sigma_z.size() == static_cast<unsigned>(nztrue_)
                        ? sigma_z.at(iz)
                        : sigma_z.at(iz + nztrue_ * ie);
            std::normal_distribution<> e{0, sigma};
            observed_events_.at(iz + rdata.nztrue * ie)
                = rdata.z.at(iz + rdata.nz * ie) + e(gen);
            observed_events_std_dev_.at(iz + rdata.nztrue * ie) = sigma;
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

void ExpData::set_observed_data(std::vector<realtype> const& observedData) {
    check_data_dimension(observedData, "observedData");

    if (observedData.size() == static_cast<unsigned>(nt()) * nytrue_)
        observed_data_ = observedData;
    else if (observedData.empty())
        observed_data_.clear();
}

void ExpData::set_observed_data(
    std::vector<realtype> const& observedData, int iy
) {
    if (observedData.size() != static_cast<unsigned>(nt()))
        throw AmiException(
            "Input observedData did not match dimensions nt (%i), was %i", nt(),
            observedData.size()
        );

    for (int it = 0; it < nt(); ++it)
        observed_data_.at(iy + it * nytrue_) = observedData.at(it);
}

bool ExpData::is_set_observed_data(int const it, int const iy) const {
    return !observed_data_.empty()
           && !isnan(observed_data_.at(it * nytrue_ + iy));
}

std::vector<realtype> const& ExpData::get_observed_data() const {
    return observed_data_;
}

realtype const* ExpData::get_observed_data_ptr(int const it) const {
    if (!observed_data_.empty())
        return &observed_data_.at(it * nytrue_);

    return nullptr;
}

void ExpData::set_observed_data_std_dev(
    std::vector<realtype> const& observedDataStdDev
) {
    check_data_dimension(observedDataStdDev, "observedDataStdDev");
    checkSigmaPositivity(observedDataStdDev, "observedDataStdDev");

    if (observedDataStdDev.size() == static_cast<unsigned>(nt()) * nytrue_)
        observed_data_std_dev_ = observedDataStdDev;
    else if (observedDataStdDev.empty())
        observed_data_std_dev_.clear();
}

void ExpData::set_observed_data_std_dev(realtype const stdDev) {
    checkSigmaPositivity(stdDev, "stdDev");
    std::ranges::fill(observed_data_std_dev_, stdDev);
}

void ExpData::set_observed_data_std_dev(
    std::vector<realtype> const& observedDataStdDev, int const iy
) {
    if (observedDataStdDev.size() != static_cast<unsigned>(nt()))
        throw AmiException(
            "Input observedDataStdDev did not match dimensions nt (%i), was %i",
            nt(), observedDataStdDev.size()
        );
    checkSigmaPositivity(observedDataStdDev, "observedDataStdDev");

    for (int it = 0; it < nt(); ++it)
        observed_data_std_dev_.at(iy + it * nytrue_)
            = observedDataStdDev.at(it);
}

void ExpData::set_observed_data_std_dev(realtype const stdDev, int const iy) {
    checkSigmaPositivity(stdDev, "stdDev");
    for (int it = 0; it < nt(); ++it)
        observed_data_std_dev_.at(iy + it * nytrue_) = stdDev;
}

bool ExpData::is_set_observed_data_std_dev(int const it, int const iy) const {
    return !observed_data_std_dev_.empty()
           && !isnan(observed_data_std_dev_.at(it * nytrue_ + iy));
}

std::vector<realtype> const& ExpData::get_observed_data_std_dev() const {
    return observed_data_std_dev_;
}

realtype const* ExpData::get_observed_data_std_dev_ptr(int const it) const {
    if (!observed_data_std_dev_.empty())
        return &observed_data_std_dev_.at(it * nytrue_);

    return nullptr;
}

void ExpData::set_observed_events(std::vector<realtype> const& observedEvents) {
    check_events_dimension(observedEvents, "observedEvents");

    if (observedEvents.size() == static_cast<unsigned>(nmaxevent_) * nztrue_)
        observed_events_ = observedEvents;
    else if (observedEvents.empty())
        observed_events_.clear();
}

void ExpData::set_observed_events(
    std::vector<realtype> const& observedEvents, int const iz
) {
    if (observedEvents.size() != static_cast<unsigned>(nmaxevent_)) {
        throw AmiException(
            "Input observedEvents did not match dimensions nmaxevent (%i), was "
            "%i",
            nmaxevent_, observedEvents.size()
        );
    }

    for (int ie = 0; ie < nmaxevent_; ++ie)
        observed_events_.at(iz + ie * nztrue_) = observedEvents.at(ie);
}

bool ExpData::is_set_observed_events(int const ie, int const iz) const {
    return !observed_events_.empty()
           && !isnan(observed_events_.at(ie * nztrue_ + iz));
}

std::vector<realtype> const& ExpData::get_observed_events() const {
    return observed_events_;
}

realtype const* ExpData::get_observed_events_ptr(int const ie) const {
    if (!observed_events_.empty())
        return &observed_events_.at(ie * nztrue_);

    return nullptr;
}

void ExpData::set_observed_events_std_dev(
    std::vector<realtype> const& observedEventsStdDev
) {
    check_events_dimension(observedEventsStdDev, "observedEventsStdDev");
    checkSigmaPositivity(observedEventsStdDev, "observedEventsStdDev");

    if (observedEventsStdDev.size() == (unsigned)nmaxevent_ * nztrue_)
        observed_events_std_dev_ = observedEventsStdDev;
    else if (observedEventsStdDev.empty())
        observed_events_std_dev_.clear();
}

void ExpData::set_observed_events_std_dev(realtype const stdDev) {
    checkSigmaPositivity(stdDev, "stdDev");
    std::ranges::fill(observed_events_std_dev_, stdDev);
}

void ExpData::set_observed_events_std_dev(
    std::vector<realtype> const& observedEventsStdDev, int const iz
) {
    if (observedEventsStdDev.size() != (unsigned)nmaxevent_)
        throw AmiException(
            "Input observedEventsStdDev did not match dimensions nmaxevent "
            "(%i), was %i",
            nmaxevent_, observedEventsStdDev.size()
        );
    checkSigmaPositivity(observedEventsStdDev, "observedEventsStdDev");

    for (int ie = 0; ie < nmaxevent_; ++ie)
        observed_events_std_dev_.at(iz + ie * nztrue_)
            = observedEventsStdDev.at(ie);
}

void ExpData::set_observed_events_std_dev(realtype const stdDev, int const iz) {
    checkSigmaPositivity(stdDev, "stdDev");

    for (int ie = 0; ie < nmaxevent_; ++ie)
        observed_events_std_dev_.at(iz + ie * nztrue_) = stdDev;
}

bool ExpData::is_set_observed_events_std_dev(int const ie, int const iz) const {
    if (!observed_events_std_dev_.empty()) // avoid out of bound memory access
        return !isnan(observed_events_std_dev_.at(ie * nztrue_ + iz));

    return false;
}

std::vector<realtype> const& ExpData::get_observed_events_std_dev() const {
    return observed_events_std_dev_;
}

realtype const* ExpData::get_observed_events_std_dev_ptr(int const ie) const {
    if (!observed_events_std_dev_.empty())
        return &observed_events_std_dev_.at(ie * nztrue_);

    return nullptr;
}

void ExpData::clear_observations() {
    std::ranges::fill(observed_data_, get_nan());
    std::ranges::fill(observed_data_std_dev_, get_nan());
    std::ranges::fill(observed_events_, get_nan());
    std::ranges::fill(observed_events_std_dev_, get_nan());
}

void ExpData::apply_dimensions() {
    apply_data_dimension();
    apply_event_dimension();
}

void ExpData::apply_data_dimension() {
    observed_data_.resize(nt() * nytrue_, get_nan());
    observed_data_std_dev_.resize(nt() * nytrue_, get_nan());
}

void ExpData::apply_event_dimension() {
    observed_events_.resize(nmaxevent_ * nztrue_, get_nan());
    observed_events_std_dev_.resize(nmaxevent_ * nztrue_, get_nan());
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

void checkSigmaPositivity(
    std::vector<realtype> const& sigmaVector, char const* vectorName
) {
    for (auto&& sigma : sigmaVector)
        checkSigmaPositivity(sigma, vectorName);
}

void checkSigmaPositivity(realtype const sigma, char const* sigmaName) {
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
    , original_parameters_(model->get_parameters())
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
    if (!edata->parameters.empty()) {
        if (edata->parameters.size() != (unsigned)model_->np())
            throw AmiException(
                "Number of parameters (%d) in model does not"
                " match ExpData (%zd).",
                model_->np(), edata->parameters.size()
            );
        model_->set_parameters(edata->parameters);
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

    model_->set_parameters(original_parameters_);
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
