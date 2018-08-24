#include "amici/edata.h"
#include "amici/rdata.h"

#include "amici/defines.h"
#include "amici/model.h"

#include <cstring>
#include <random>
#include <utility>
#include <algorithm>

namespace amici {

ExpData::ExpData() : nytrue(0), nztrue(0), nmaxevent(0) {}

ExpData::ExpData(int nytrue, int nztrue, int nmaxevent)
    : nytrue(nytrue), nztrue(nztrue), nmaxevent(nmaxevent)
{
}

ExpData::ExpData(int nytrue, int nztrue, int nmaxevent,
                 std::vector<realtype> ts)
    : nytrue(nytrue), nztrue(nztrue), nmaxevent(nmaxevent)
{
    setTimepoints(ts);
}

ExpData::ExpData(int nytrue, int nztrue, int nmaxevent,
                 std::vector<realtype> ts,
                 std::vector<realtype> observedData,
                 std::vector<realtype> observedDataStdDev,
                 std::vector<realtype> observedEvents,
                 std::vector<realtype> observedEventsStdDev)
    : nytrue(nytrue), nztrue(nztrue), nmaxevent(nmaxevent), ts(std::move(ts))
{
    setObservedData(observedData);
    setObservedDataStdDev(observedDataStdDev);
    setObservedEvents(observedEvents);
    setObservedEventsStdDev(observedEventsStdDev);
}

ExpData::ExpData(Model const& model)
    : ExpData(model.nytrue, model.nztrue, model.nMaxEvent(), model.getTimepoints())
{
    fixedParameters = std::move(model.getFixedParameters());
}

ExpData::ExpData(const ExpData &other)
    : ExpData(other.nytrue, other.nztrue, other.nmaxevent,
              other.ts, other.observedData, other.observedDataStdDev, other.observedEvents, other.observedEventsStdDev)
{
    fixedParameters = std::move(other.fixedParameters);
    fixedParametersPreequilibration = std::move(other.fixedParametersPreequilibration);
}
    
ExpData::ExpData(ReturnData const& rdata, realtype sigma_y, realtype sigma_z)
    : ExpData(rdata, std::vector<realtype>(rdata.nytrue*rdata.nt, sigma_y), std::vector<realtype>(rdata.nztrue*rdata.nmaxevent, sigma_z))
{
}
    
ExpData::ExpData(ReturnData const& rdata, std::vector<realtype> sigma_y, std::vector<realtype> sigma_z)
    : ExpData(rdata.nytrue, rdata.nztrue, rdata.nmaxevent, rdata.ts)
{
    if (sigma_y.size() != (unsigned) nytrue && sigma_y.size() != (unsigned) nytrue*nt())
        throw AmiException("Dimension of sigma_y must be %d or %d, was %d", nytrue, nytrue*nt(), sigma_y.size());
    
    if (sigma_z.size() != (unsigned) nztrue && sigma_z.size() != (unsigned) nztrue*nmaxevent)
        throw AmiException("Dimension of sigma_z must be %d or %d, was %d", nztrue, nztrue*nmaxevent, sigma_z.size());
    
    std::random_device rd{};
    std::mt19937 gen{rd()};
    
    realtype sigma;
    
    checkSigmaPositivity(sigma_y, "sigma_y");
    checkSigmaPositivity(sigma_z, "sigma_z");
    
    for (int iy = 0; iy < nytrue; ++iy) {
        for (int it = 0; it < nt(); ++it) {
            sigma = sigma_y.size() == (unsigned) nytrue ? sigma_y.at(iy) : sigma_y.at(iy + nytrue * it);
            std::normal_distribution<> e{0, sigma};
            observedData.at(iy + nytrue * it) = rdata.y.at(iy + rdata.ny * it) + e(gen);
            observedDataStdDev.at(iy + nytrue * it) = sigma;
        }
    }
    
    for (int iz = 0; iz < nztrue; ++iz) {
        for (int ie = 0; ie < nmaxevent; ++ie) {
            sigma = sigma_z.size() == (unsigned) nztrue ? sigma_z.at(iz) : sigma_z.at(iz + nztrue * ie);
            std::normal_distribution<> e{0, sigma};
            observedEvents.at(iz + rdata.nztrue * ie) = rdata.z.at(iz + rdata.nz * ie) + e(gen);
            observedDataStdDev.at(iz + rdata.nztrue * ie) = sigma;
        }
        }
}
    
void ExpData::setTimepoints(const std::vector<realtype> &ts) {
    if (!std::is_sorted(ts.begin(), ts.end()))
        throw AmiException("Encountered non-monotonic timepoints, please order timepoints such that they are monotonically increasing!");

    this->ts = std::move(ts);
    observedData.resize(nt()*nytrue, getNaN());
    observedDataStdDev.resize(nt()*nytrue, getNaN());
}
    
std::vector<realtype> const& ExpData::getTimepoints() const {
    return ts;
}

const int ExpData::nt() const {
    return ts.size();
}
    
realtype ExpData::getTimepoint(int it) const {
    return ts.at(it);
}
    
void ExpData::setObservedData(const std::vector<realtype> &observedData) {
    checkDataDimension(observedData, "observedData");
        
    if (observedData.size() == (unsigned) nt()*nytrue)
        this->observedData = std::move(observedData);
    else if (observedData.empty())
        this->observedData.clear();
}
    
void ExpData::setObservedData(const std::vector<realtype> &observedData, int iy) {
    if (observedData.size() != (unsigned) nt())
        throw AmiException("Input observedData did not match dimensions nt (%i), was %i", nt(), observedData.size());
    
    for (int it = 0; it < nt(); ++it)
        this->observedData.at(iy + it*nytrue) = observedData.at(it);
}
    
bool ExpData::isSetObservedData(int it, int iy) const {
    return !observedData.empty() && !isNaN(observedData.at(it * nytrue + iy));
}
    
std::vector<realtype> const& ExpData::getObservedData() const {
    return observedData;
}
    
const realtype *ExpData::getObservedDataPtr(int it) const {
    if (!observedData.empty())
        return &observedData.at(it*nytrue);
    
    return nullptr;
}

void ExpData::setObservedDataStdDev(const std::vector<realtype> &observedDataStdDev) {
    checkDataDimension(observedDataStdDev, "observedDataStdDev");
    checkSigmaPositivity(observedDataStdDev, "observedDataStdDev");
    
    if (observedDataStdDev.size() == (unsigned) nt()*nytrue)
        this->observedDataStdDev = std::move(observedDataStdDev);
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
        this->observedDataStdDev.at(iy + it*nytrue) = observedDataStdDev.at(it);
}
    
void ExpData::setObservedDataStdDev(const realtype stdDev, int iy) {
    checkSigmaPositivity(stdDev, "stdDev");
    for (int it = 0; it < nt(); ++it)
        observedDataStdDev.at(iy + it*nytrue) = stdDev;
}
    
bool ExpData::isSetObservedDataStdDev(int it, int iy) const {
    return !observedDataStdDev.empty() && !isNaN(observedDataStdDev.at(it * nytrue + iy));
}
    
std::vector<realtype> const& ExpData::getObservedDataStdDev() const {
    return observedDataStdDev;
}

const realtype *ExpData::getObservedDataStdDevPtr(int it) const {
    if (!observedDataStdDev.empty())
        return &observedDataStdDev.at(it*nytrue);
    
    return nullptr;
}
    
void ExpData::setObservedEvents(const std::vector<realtype> &observedEvents) {
    checkEventsDimension(observedEvents, "observedEvents");
    
    if (observedEvents.size() == (unsigned) nmaxevent*nztrue)
        this->observedEvents = std::move(observedEvents);
    else if (observedEvents.empty())
        this->observedEvents.clear();
}
    
void ExpData::setObservedEvents(const std::vector<realtype> &observedEvents, int iz) {
    if (observedEvents.size() != (unsigned) nmaxevent) {
        throw AmiException("Input observedEvents did not match dimensions nmaxevent (%i), was %i", nmaxevent, observedEvents.size());
    }
    
    for (int ie = 0; ie < nmaxevent; ++ie)
        this->observedEvents.at(iz + ie*nztrue) = observedEvents.at(ie);
}
    
bool ExpData::isSetObservedEvents(int ie, int iz) const {
    return !observedEvents.size() && !isNaN(observedEvents.at(ie * nztrue + iz));
}

std::vector<realtype> const& ExpData::getObservedEvents() const {
    return observedEvents;
}

const realtype *ExpData::getObservedEventsPtr(int ie) const {
    if (!observedEvents.empty())
        return &observedEvents.at(ie*nztrue);
    
    return nullptr;
}
    
void ExpData::setObservedEventsStdDev(const std::vector<realtype> &observedEventsStdDev) {
    checkEventsDimension(observedEventsStdDev, "observedEventsStdDev");
    checkSigmaPositivity(observedEventsStdDev, "observedEventsStdDev");
        
    if (observedEventsStdDev.size() == (unsigned) nmaxevent*nztrue)
        this->observedEventsStdDev = std::move(observedEventsStdDev);
    else if (observedEventsStdDev.empty())
        this->observedEventsStdDev.clear();
}
    
void ExpData::setObservedEventsStdDev(const realtype stdDev) {
    checkSigmaPositivity(stdDev, "stdDev");
    std::fill(observedEventsStdDev.begin() ,observedEventsStdDev.end(), stdDev);
}

void ExpData::setObservedEventsStdDev(const std::vector<realtype> &observedEventsStdDev, int iz) {
    if (observedEventsStdDev.size() != (unsigned) nmaxevent)
        throw AmiException("Input observedEventsStdDev did not match dimensions nmaxevent (%i), was %i", nmaxevent, observedEventsStdDev.size());
    checkSigmaPositivity(observedEventsStdDev, "observedEventsStdDev");
        
    for (int ie = 0; ie < nmaxevent; ++ie)
        this->observedEventsStdDev.at(iz + ie*nztrue) = observedEventsStdDev.at(ie);
}

void ExpData::setObservedEventsStdDev(const realtype stdDev, int iz) {
    checkSigmaPositivity(stdDev, "stdDev");
    
    for (int ie = 0; ie < nmaxevent; ++ie)
        observedEventsStdDev.at(iz + ie*nztrue) = stdDev;
}
    
bool ExpData::isSetObservedEventsStdDev(int ie, int iz) const {
    if (!observedEventsStdDev.empty()) // avoid out of bound memory access
        return !isNaN(observedEventsStdDev.at(ie * nztrue + iz));
  
    return false;
}

std::vector<realtype> const& ExpData::getObservedEventsStdDev() const {
    return observedEventsStdDev;
}

const realtype *ExpData::getObservedEventsStdDevPtr(int ie) const {
    if (!observedEventsStdDev.empty())
        return &observedEventsStdDev.at(ie*nztrue);
   
    return nullptr;
}
    
void ExpData::checkDataDimension(std::vector<realtype> input, const char *fieldname) const {
    if (input.size() != (unsigned) nt()*nytrue && !input.empty())
        throw AmiException("Input %s did not match dimensions nt (%i) x nytrue (%i), was %i", fieldname, nt(), nytrue, input.size());
}
    
void ExpData::checkEventsDimension(std::vector<realtype> input, const char *fieldname) const {
    if (input.size() != (unsigned) nmaxevent*nztrue && !input.empty())
        throw AmiException("Input %s did not match dimensions nt (%i) x nytrue (%i), was %i", fieldname, nmaxevent, nztrue, input.size());
}

void ExpData::checkSigmaPositivity(std::vector<realtype> sigmaVector, const char *vectorName) const {
    for (auto&& sigma : sigmaVector)
        checkSigmaPositivity(sigma, vectorName);
}

void ExpData::checkSigmaPositivity(realtype sigma, const char *sigmaName) const {
    if (sigma <= 0.0)
        throw AmiException("Encountered sigma <= 0 in %s! value: %f", sigmaName, sigma);
}

} // namespace amici
