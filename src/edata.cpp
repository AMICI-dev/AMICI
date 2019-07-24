#include "amici/edata.h"
#include "amici/rdata.h"
#include "amici/symbolic_functions.h" // getNaN
#include "amici/defines.h"
#include "amici/model.h"

#include <cstring>
#include <random>
#include <utility>
#include <algorithm>

namespace amici {

ExpData::ExpData(int nytrue, int nztrue, int nmaxevent)
    : nytrue_(nytrue), nztrue_(nztrue), nmaxevent_(nmaxevent)
{
    applyDimensions();
}

ExpData::ExpData(int nytrue, int nztrue, int nmaxevent,
                 std::vector<realtype> ts_)
    : nytrue_(nytrue), nztrue_(nztrue), nmaxevent_(nmaxevent), ts(std::move(ts_))
{
    applyDimensions();
}

ExpData::ExpData(int nytrue, int nztrue, int nmaxevent,
                 std::vector<realtype> ts_,
                 std::vector<realtype> fixedParameters_
                 )
    : fixedParameters(std::move(fixedParameters_)), nytrue_(nytrue), nztrue_(nztrue), nmaxevent_(nmaxevent), ts(std::move(ts_))
{
    applyDimensions();
}

ExpData::ExpData(int nytrue, int nztrue, int nmaxevent,
                 std::vector<realtype> ts_,
                 std::vector<realtype> const& observedData,
                 std::vector<realtype> const& observedDataStdDev,
                 std::vector<realtype> const& observedEvents,
                 std::vector<realtype> const& observedEventsStdDev)
    : nytrue_(nytrue), nztrue_(nztrue), nmaxevent_(nmaxevent), ts(std::move(ts_))
{
    applyDimensions();
    setObservedData(observedData);
    setObservedDataStdDev(observedDataStdDev);
    setObservedEvents(observedEvents);
    setObservedEventsStdDev(observedEventsStdDev);
}

ExpData::ExpData(Model const& model)
    : ExpData(model.nytrue, model.nztrue, model.nMaxEvent(),
              model.getTimepoints(), model.getFixedParameters()) {}

ExpData::ExpData(ReturnData const& rdata, realtype sigma_y, realtype sigma_z)
    : ExpData(rdata, std::vector<realtype>(rdata.nytrue*rdata.nt, sigma_y), std::vector<realtype>(rdata.nztrue*rdata.nmaxevent, sigma_z)) {}

ExpData::ExpData(ReturnData const& rdata, std::vector<realtype> sigma_y,
                 std::vector<realtype> sigma_z)
    : ExpData(rdata.nytrue, rdata.nztrue, rdata.nmaxevent, rdata.ts)
{
    if (sigma_y.size() != (unsigned) nytrue_ && sigma_y.size() != (unsigned) nytrue_*nt())
        throw AmiException("Dimension of sigma_y must be %d or %d, was %d", nytrue_, nytrue_*nt(), sigma_y.size());

    if (sigma_z.size() != (unsigned) nztrue_ && sigma_z.size() != (unsigned) nztrue_*nmaxevent_)
        throw AmiException("Dimension of sigma_z must be %d or %d, was %d", nztrue_, nztrue_*nmaxevent_, sigma_z.size());

    std::random_device rd{};
    std::mt19937 gen{rd()};

    realtype sigma;

    checkSigmaPositivity(sigma_y, "sigma_y");
    checkSigmaPositivity(sigma_z, "sigma_z");

    for (int iy = 0; iy < nytrue_; ++iy) {
        for (int it = 0; it < nt(); ++it) {
            sigma = sigma_y.size() == (unsigned) nytrue_ ? sigma_y.at(iy) : sigma_y.at(iy + nytrue_ * it);
            std::normal_distribution<> e{0, sigma};
            observedData.at(iy + nytrue_ * it) = rdata.y.at(iy + rdata.ny * it) + e(gen);
            observedDataStdDev.at(iy + nytrue_ * it) = sigma;
        }
    }

    for (int iz = 0; iz < nztrue_; ++iz) {
        for (int ie = 0; ie < nmaxevent_; ++ie) {
            sigma = sigma_z.size() == (unsigned) nztrue_ ? sigma_z.at(iz) : sigma_z.at(iz + nztrue_ * ie);
            std::normal_distribution<> e{0, sigma};
            observedEvents.at(iz + rdata.nztrue * ie) = rdata.z.at(iz + rdata.nz * ie) + e(gen);
            observedDataStdDev.at(iz + rdata.nztrue * ie) = sigma;
        }
        }
}

void ExpData::setTimepoints(const std::vector<realtype> &ts) {
    if (!std::is_sorted(ts.begin(), ts.end()))
        throw AmiException("Encountered non-monotonic timepoints, please order timepoints such that they are monotonically increasing!");
    this->ts = ts;
    applyDataDimension();
}

std::vector<realtype> const& ExpData::getTimepoints() const {
    return ts;
}

int ExpData::nt() const {
    return ts.size();
}

realtype ExpData::getTimepoint(int it) const {
    return ts.at(it);
}

void ExpData::setObservedData(const std::vector<realtype> &observedData) {
    checkDataDimension(observedData, "observedData");

    if (observedData.size() == (unsigned) nt()*nytrue_)
        this->observedData = observedData;
    else if (observedData.empty())
        this->observedData.clear();
}

void ExpData::setObservedData(const std::vector<realtype> &observedData, int iy) {
    if (observedData.size() != (unsigned) nt())
        throw AmiException("Input observedData did not match dimensions nt (%i), was %i", nt(), observedData.size());

    for (int it = 0; it < nt(); ++it)
        this->observedData.at(iy + it*nytrue_) = observedData.at(it);
}

bool ExpData::isSetObservedData(int it, int iy) const {
    return !observedData.empty() && !isNaN(observedData.at(it * nytrue_ + iy));
}

std::vector<realtype> const& ExpData::getObservedData() const {
    return observedData;
}

const realtype *ExpData::getObservedDataPtr(int it) const {
    if (!observedData.empty())
        return &observedData.at(it*nytrue_);

    return nullptr;
}

void ExpData::setObservedDataStdDev(const std::vector<realtype> &observedDataStdDev) {
    checkDataDimension(observedDataStdDev, "observedDataStdDev");
    checkSigmaPositivity(observedDataStdDev, "observedDataStdDev");

    if (observedDataStdDev.size() == (unsigned) nt()*nytrue_)
        this->observedDataStdDev = observedDataStdDev;
    else if (observedDataStdDev.empty())
        this->observedDataStdDev.clear();
}

void ExpData::setObservedDataStdDev(const realtype stdDev) {
    checkSigmaPositivity(stdDev, "stdDev");
    std::fill(observedDataStdDev.begin() ,observedDataStdDev.end(), stdDev);
}

void ExpData::setObservedDataStdDev(const std::vector<realtype> &observedDataStdDev, int iy) {
    if (observedDataStdDev.size() != (unsigned) nt())
        throw AmiException("Input observedDataStdDev did not match dimensions nt (%i), was %i", nt(), observedDataStdDev.size());
    checkSigmaPositivity(observedDataStdDev, "observedDataStdDev");

    for (int it = 0; it < nt(); ++it)
        this->observedDataStdDev.at(iy + it*nytrue_) = observedDataStdDev.at(it);
}

void ExpData::setObservedDataStdDev(const realtype stdDev, int iy) {
    checkSigmaPositivity(stdDev, "stdDev");
    for (int it = 0; it < nt(); ++it)
        observedDataStdDev.at(iy + it*nytrue_) = stdDev;
}

bool ExpData::isSetObservedDataStdDev(int it, int iy) const {
    return !observedDataStdDev.empty() && !isNaN(observedDataStdDev.at(it * nytrue_ + iy));
}

std::vector<realtype> const& ExpData::getObservedDataStdDev() const {
    return observedDataStdDev;
}

const realtype *ExpData::getObservedDataStdDevPtr(int it) const {
    if (!observedDataStdDev.empty())
        return &observedDataStdDev.at(it*nytrue_);

    return nullptr;
}

void ExpData::setObservedEvents(const std::vector<realtype> &observedEvents) {
    checkEventsDimension(observedEvents, "observedEvents");

    if (observedEvents.size() == (unsigned) nmaxevent_*nztrue_)
        this->observedEvents = observedEvents;
    else if (observedEvents.empty())
        this->observedEvents.clear();
}

void ExpData::setObservedEvents(const std::vector<realtype> &observedEvents, int iz) {
    if (observedEvents.size() != (unsigned) nmaxevent_) {
        throw AmiException("Input observedEvents did not match dimensions nmaxevent (%i), was %i", nmaxevent_, observedEvents.size());
    }

    for (int ie = 0; ie < nmaxevent_; ++ie)
        this->observedEvents.at(iz + ie*nztrue_) = observedEvents.at(ie);
}

bool ExpData::isSetObservedEvents(int ie, int iz) const {
    return !observedEvents.empty() && !isNaN(observedEvents.at(ie * nztrue_ + iz));
}

std::vector<realtype> const& ExpData::getObservedEvents() const {
    return observedEvents;
}

const realtype *ExpData::getObservedEventsPtr(int ie) const {
    if (!observedEvents.empty())
        return &observedEvents.at(ie*nztrue_);

    return nullptr;
}

void ExpData::setObservedEventsStdDev(const std::vector<realtype> &observedEventsStdDev) {
    checkEventsDimension(observedEventsStdDev, "observedEventsStdDev");
    checkSigmaPositivity(observedEventsStdDev, "observedEventsStdDev");

    if (observedEventsStdDev.size() == (unsigned) nmaxevent_*nztrue_)
        this->observedEventsStdDev = observedEventsStdDev;
    else if (observedEventsStdDev.empty())
        this->observedEventsStdDev.clear();
}

void ExpData::setObservedEventsStdDev(const realtype stdDev) {
    checkSigmaPositivity(stdDev, "stdDev");
    std::fill(observedEventsStdDev.begin() ,observedEventsStdDev.end(), stdDev);
}

void ExpData::setObservedEventsStdDev(const std::vector<realtype> &observedEventsStdDev, int iz) {
    if (observedEventsStdDev.size() != (unsigned) nmaxevent_)
        throw AmiException("Input observedEventsStdDev did not match dimensions nmaxevent (%i), was %i", nmaxevent_, observedEventsStdDev.size());
    checkSigmaPositivity(observedEventsStdDev, "observedEventsStdDev");

    for (int ie = 0; ie < nmaxevent_; ++ie)
        this->observedEventsStdDev.at(iz + ie*nztrue_) = observedEventsStdDev.at(ie);
}

void ExpData::setObservedEventsStdDev(const realtype stdDev, int iz) {
    checkSigmaPositivity(stdDev, "stdDev");

    for (int ie = 0; ie < nmaxevent_; ++ie)
        observedEventsStdDev.at(iz + ie*nztrue_) = stdDev;
}

bool ExpData::isSetObservedEventsStdDev(int ie, int iz) const {
    if (!observedEventsStdDev.empty()) // avoid out of bound memory access
        return !isNaN(observedEventsStdDev.at(ie * nztrue_ + iz));

    return false;
}

std::vector<realtype> const& ExpData::getObservedEventsStdDev() const {
    return observedEventsStdDev;
}

const realtype *ExpData::getObservedEventsStdDevPtr(int ie) const {
    if (!observedEventsStdDev.empty())
        return &observedEventsStdDev.at(ie*nztrue_);

    return nullptr;
}

void ExpData::applyDimensions() {
    applyDataDimension();
    applyEventDimension();
}

void ExpData::applyDataDimension() {
    observedData.resize(nt()*nytrue_, getNaN());
    observedDataStdDev.resize(nt()*nytrue_, getNaN());
}

void ExpData::applyEventDimension() {
    observedEvents.resize(nmaxevent_*nztrue_, getNaN());
    observedEventsStdDev.resize(nmaxevent_*nztrue_, getNaN());
}

void ExpData::checkDataDimension(std::vector<realtype> const& input, const char *fieldname) const {
    if (input.size() != (unsigned) nt()*nytrue_ && !input.empty())
        throw AmiException("Input %s did not match dimensions nt (%i) x nytrue (%i), was %i", fieldname, nt(), nytrue_, input.size());
}

void ExpData::checkEventsDimension(std::vector<realtype> const& input, const char *fieldname) const {
    if (input.size() != (unsigned) nmaxevent_*nztrue_ && !input.empty())
        throw AmiException("Input %s did not match dimensions nt (%i) x nytrue (%i), was %i", fieldname, nmaxevent_, nztrue_, input.size());
}

void checkSigmaPositivity(std::vector<realtype> const& sigmaVector, const char *vectorName) {
    for (auto&& sigma : sigmaVector)
        checkSigmaPositivity(sigma, vectorName);
}

void checkSigmaPositivity(const realtype sigma, const char *sigmaName) {
    if (sigma <= 0.0)
        throw AmiException("Encountered sigma <= 0 in %s! value: %f", sigmaName, sigma);
}

int ExpData::nytrue() const
{
    return nytrue_;
}

int ExpData::nztrue() const
{
    return nztrue_;
}

int ExpData::nmaxevent() const
{
    return nmaxevent_;
}

ConditionContext::ConditionContext(Model *model, const ExpData *edata)
    : model(model),
      originalx0(model->getInitialStates()),
      originalsx0(model->getInitialStateSensitivities()),
      originalParameters(model->getParameters()),
      originalFixedParameters(model->getFixedParameters()),
      originalTimepoints(model->getTimepoints()),
      originalParameterList(model->getParameterList()),
      originalScaling(model->getParameterScale())
{
    applyCondition(edata);
}

ConditionContext::~ConditionContext()
{
    restore();
}

void ConditionContext::applyCondition(const ExpData *edata)
{
    if(!edata)
        return;

    // this needs to go first, otherwise nplist will not have the right
    // dimension for all other fields that depend on Model::nplist
    if(!edata->plist.empty())
        model->setParameterList(edata->plist);

    // this needs to go second as setParameterScale will reset sx0
    if(!edata->pscale.empty()) {
        if(edata->pscale.size() != (unsigned) model->np())
            throw AmiException("Number of parameters (%d) in model does not"
                               " match ExpData (%zd).",
                               model->np(), edata->pscale.size());
        model->setParameterScale(edata->pscale);

    }

    if(!edata->x0.empty()) {
        if(edata->x0.size() != (unsigned) model->nx_rdata)
            throw AmiException("Number of initial conditions (%d) in model does"
                               " not match ExpData (%zd).",
                               model->nx_rdata, edata->x0.size());
        model->setInitialStates(edata->x0);
    }

    if(!edata->sx0.empty()) {
        if(edata->sx0.size() != (unsigned) model->nx_rdata * model->nplist())
            throw AmiException("Number of initial conditions sensitivities (%d)"
                               " in model does not match ExpData (%zd).",
                               model->nx_rdata * model->nplist(),
                               edata->sx0.size());
        model->setInitialStateSensitivities(edata->sx0);
    }

    if(!edata->parameters.empty()) {
        if(edata->parameters.size() != (unsigned) model->np())
            throw AmiException("Number of parameters (%d) in model does not"
                               " match ExpData (%zd).",
                               model->np(), edata->parameters.size());
        model->setParameters(edata->parameters);
    }

    if(!edata->fixedParameters.empty()) {
        // fixed parameter in model are superseded by those provided in edata
        if(edata->fixedParameters.size() != (unsigned) model->nk())
            throw AmiException("Number of fixed parameters (%d) in model does"
                               " not match ExpData (%zd).",
                               model->nk(), edata->fixedParameters.size());
        model->setFixedParameters(edata->fixedParameters);
    }

    if(edata->nt()) {
        // fixed parameter in model are superseded by those provided in edata
        model->setTimepoints(edata->getTimepoints());
    }
}

void ConditionContext::restore()
{
    // parameter list has to be set before initial state sensitivities
    model->setParameterList(originalParameterList);
    // parameter scale has to be set before initial state sensitivities
    model->setParameterScale(originalScaling);
    model->setInitialStates(originalx0);
    model->setUnscaledInitialStateSensitivities(originalsx0);
    model->setParameters(originalParameters);
    model->setFixedParameters(originalFixedParameters);
    model->setTimepoints(originalTimepoints);

}


} // namespace amici
