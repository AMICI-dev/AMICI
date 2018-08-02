#include "amici/edata.h"
#include "amici/rdata.h"

#include "amici/defines.h"
#include "amici/model.h"

#include <cstring>
#include <random>
#include <utility>

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
                 std::vector<realtype> my,
                 std::vector<realtype> sigmay,
                 std::vector<realtype> mz,
                 std::vector<realtype> sigmaz)
    : nytrue(nytrue), nztrue(nztrue), nmaxevent(nmaxevent), ts(std::move(ts)), observedData(std::move(my)), observedDataStdDev(std::move(sigmay)), observedEvents(std::move(mz)), observedEventsStdDev(std::move(sigmaz))
{
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
    
    for (int iy = 0; iy < nytrue; ++iy) {
        for (int it = 0; it < nt(); ++it) {
            if (sigma_y.size() == (unsigned) nytrue) {
                if (sigma_y.at(iy) < 0.0)
                    throw AmiException("All sigma_y must positive! sigma at index %d was %f", iy, sigma_y[iy]);
                sigma = sigma_y.at(iy);
            } else {
                if (sigma_y.at(iy + nytrue * it) < 0.0)
                    throw AmiException("All sigma_y must positive! sigma at index %d was %f", iy + nytrue * it, sigma_y[iy + nytrue * it]);
                sigma = sigma_y.at(iy + nytrue * it);
            }
            std::normal_distribution<> e{0, sigma};
            observedData.at(iy + nytrue * it) = rdata.y.at(iy + rdata.ny * it) + e(gen);
            observedDataStdDev.at(iy + nytrue * it) = sigma;
        }
    }
    
    for (int iz = 0; iz < nztrue; ++iz) {
        for (int ie = 0; ie < nmaxevent; ++ie) {
            if (sigma_z.size() == (unsigned) nztrue) {
                if (sigma_z.at(iz) < 0.0)
                    throw AmiException("All sigma_z must positive! sigma at index %d was %f", iz, sigma_z.at(iz));
                sigma = sigma_z.at(iz);
            } else {
                if (sigma_z.at(iz + rdata.nztrue * ie) < 0.0)
                    throw AmiException("All sigma_z must positive! sigma at index %d was %f", iz + rdata.nztrue * ie, sigma_z.at(iz + rdata.nztrue * ie));
                sigma = sigma_z.at(iz + rdata.nztrue * ie);
            }
            std::normal_distribution<> e{0, sigma};
            observedEvents.at(iz + rdata.nztrue * ie) = rdata.z.at(iz + rdata.nz * ie) + e(gen);
            observedDataStdDev.at(iz + rdata.nztrue * ie) = sigma;
        }
        }
}
    
void ExpData::setTimepoints(const std::vector<realtype> &ts) {
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
    if (observedData.size() != (unsigned) nt()*nytrue && !observedData.empty())
        throw AmiException("Input observedData did not match dimensions nmaxevent (%i) x nztrue (%i), was %i", nt(), nytrue, observedData.size());
        
    if (observedData.size() == (unsigned) nt()*nytrue)
        this->observedData = std::move(observedData);
    else if (observedData.empty())
        this->observedData.clear();
}
    
void ExpData::setObservedData(const std::vector<realtype> &observedData, int iy) {
    if (observedData.size() != (unsigned) nt())
        throw AmiException("Input observedData did not match dimensions nt (%i), was %i", nt(), observedData.size());
    
    for (int it = 0; it < nt(); ++it)
        this->observedData.at(iy + it*nytrue) = observedData.at(iy);
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
    if (observedDataStdDev.size() != (unsigned) nt()*nytrue && !observedDataStdDev.empty())
        throw AmiException("Input observedDataStdDev did not match dimensions nmaxevent (%i) x nztrue (%i), was %i", nt(), nytrue, observedDataStdDev.size());
    
    if (observedDataStdDev.size() == (unsigned) nt()*nytrue)
        this->observedDataStdDev = std::move(observedDataStdDev);
    else if (observedDataStdDev.empty())
        this->observedDataStdDev.clear();
}
    
void ExpData::setObservedDataStdDev(const realtype stdDev) {
    std::fill(observedDataStdDev.begin() ,observedDataStdDev.end(), stdDev);
}
    
void ExpData::setObservedDataStdDev(const std::vector<realtype> &observedDataStdDev, int iy) {
    if (observedDataStdDev.size() != (unsigned) nt())
        throw AmiException("Input observedDataStdDev did not match dimensions nt (%i), was %i", nt(), observedDataStdDev.size());
    
    for (int it = 0; it < nt(); ++it)
        this->observedDataStdDev.at(iy + it*nytrue) = observedDataStdDev.at(iy);
}
    
void ExpData::setObservedDataStdDev(const realtype stdDev, int iy) {
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
    if (observedEvents.size() != (unsigned) nmaxevent*nztrue && !observedEvents.empty())
        throw AmiException("Input observedEvents did not match dimensions nmaxevent (%i) x nztrue (%i), was %i", nmaxevent, nztrue, observedEvents.size());
    
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
        this->observedEvents.at(iz + ie*nztrue) = observedEvents.at(iz);
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
    if (observedEventsStdDev.size() != (unsigned) nmaxevent*nztrue && !observedEventsStdDev.empty())
        throw AmiException("Input observedEventsStdDev did not match dimensions nmaxevent (%i) x nztrue (%i), was %i", nmaxevent, nztrue, observedEventsStdDev.size());
        
    if (observedEventsStdDev.size() == (unsigned) nmaxevent*nztrue)
        this->observedEventsStdDev = std::move(observedEventsStdDev);
    else if (observedEventsStdDev.empty())
        this->observedEventsStdDev.clear();
}
    
void ExpData::setObservedEventsStdDev(const realtype stdDev) {
    std::fill(observedEventsStdDev.begin() ,observedEventsStdDev.end(), stdDev);
}

void ExpData::setObservedEventsStdDev(const std::vector<realtype> &observedEventsStdDev, int iz) {
    if (observedEventsStdDev.size() != (unsigned) nmaxevent)
        throw AmiException("Input observedEventsStdDev did not match dimensions nmaxevent (%i), was %i", nmaxevent, observedEventsStdDev.size());
        
    for (int ie = 0; ie < nmaxevent; ++ie)
        this->observedEventsStdDev.at(iz + ie*nztrue) = observedEventsStdDev.at(iz);
}

void ExpData::setObservedEventsStdDev(const realtype stdDev, int iz) {
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

} // namespace amici
