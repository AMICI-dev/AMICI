#include "amici/edata.h"
#include "amici/rdata.h"

#include "amici/defines.h"
#include "amici/model.h"

#include <cstring>
#include <random>
#include <utility>

namespace amici {

ExpData::ExpData() : nytrue(0), nztrue(0), nmaxevent(0) {}

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
    fixedParameters = model.getFixedParameters();
}

ExpData::ExpData(const ExpData &other)
    : ExpData(other.nytrue, other.nztrue, other.nmaxevent,
              other.ts, other.observedData, other.observedDataStdDev, other.observedEvents, other.observedEventsStdDev)
{
    fixedParameters = other.fixedParameters;
    fixedParametersPreequilibration = other.fixedParametersPreequilibration;
}
    
ExpData::ExpData(ReturnData const& rdata, realtype sigma_y, realtype sigma_z)
    : ExpData(rdata, std::vector<realtype>(rdata.nytrue, sigma_y), std::vector<realtype>(rdata.nztrue, sigma_z))
{
}
    
ExpData::ExpData(ReturnData const& rdata, std::vector<realtype> sigma_y, std::vector<realtype> sigma_z)
    : ExpData(rdata.nytrue, rdata.nztrue, rdata.nmaxevent, rdata.ts)
{
    if (sigma_y.size() != (unsigned) nytrue)
        throw AmiException("Dimension of sigma_y must be %d, was %d", nytrue, sigma_y.size());
        
    if (sigma_z.size() != (unsigned) nztrue)
        throw AmiException("Dimension of sigma_z must be %d, was %d", nztrue, sigma_z.size());
    
    std::random_device rd{};
    std::mt19937 gen{rd()};
    
    for (int iy = 0; iy < nytrue; ++iy) {
        if (sigma_y[iy] < 0.0)
            throw AmiException("All sigma_y must positive! sigma at index %d was %f", iy, sigma_y[iy]);
        std::normal_distribution<> e{0, sigma_y[iy]};
        for (int it = 0; it < nt(); ++it) {
            observedData.at(iy + rdata.nytrue * it) = rdata.y.at(iy + rdata.ny * it) + e(gen);
            observedDataStdDev.at(iy + rdata.ny * it) = sigma_y[iy];
        }
    }
    
    for (int iz = 0; iz < nztrue; ++iz) {
        if (sigma_z[iz] < 0.0)
            throw AmiException("All sigma_z must positive! sigma at index %d was %f", iz, sigma_z[iz]);
        std::normal_distribution<> e{0, sigma_z[iz]};
        for (int ie = 0; ie < nmaxevent; ++ie) {
            observedEvents.at(iz + rdata.nztrue * ie) = rdata.z.at(iz + rdata.nz * ie) + e(gen);
            observedDataStdDev.at(iz + rdata.nztrue * ie) = sigma_z[iz];
        }
    }
}

void ExpData::setTimepoints(const std::vector<realtype> &ts) {
    this->ts = ts;
    observedData.clear();
    observedData.resize(nt()*nytrue, getNaN());
    observedDataStdDev.clear();
    observedDataStdDev.resize(nt()*nytrue, getNaN());
}
    
std::vector<realtype> ExpData::getTimepoints() const {
    return ts;
}

const int ExpData::nt() const {
    return ts.size();
}
    
realtype ExpData::getTimepoint(int it) const {
    return ts.at(it);
}
    
void ExpData::setObservedData(const std::vector<realtype> &observedData) {
    if (observedData.size() == (unsigned) nt()*nytrue) {
        this->observedData = observedData;
    } else {
        if (observedData.size() == 0) {
            this->observedData.clear();
        } else {
            throw AmiException("Input observedData did not match dimensions nt (%i) x nytrue (%i), was %i", nt(), nytrue, observedData.size());
        }
    }
}
    
void ExpData::setObservedData(const std::vector<realtype> &observedData, int iy) {
    if (observedData.size() == (unsigned) nt()) {
        for (int it = 0; it < nt(); ++it)
            this->observedData.at(iy + it*nytrue) = observedData.at(iy);
    } else {
        throw AmiException("Input observedData did not match dimensions nt (%i), was %i", nt(), observedData.size());
    }
}
    
bool ExpData::isSetObservedData(int it, int iy) const {
    if (!observedData.empty()) {
        return !isNaN(observedData.at(it * nytrue + iy));
    } else {
        return false;
    }
}
    
std::vector<realtype> ExpData::getObservedData() const {
    return observedData;
}
    
const realtype *ExpData::getObservedDataPtr(int it) const {
    if (!observedData.empty())
        return &observedData.at(it*nytrue);
    else
        return nullptr;
}

void ExpData::setObservedDataStdDev(const std::vector<realtype> &observedDataStdDev) {
    if (observedDataStdDev.size() == (unsigned) nt()*nytrue) {
        this->observedDataStdDev = observedDataStdDev;
    } else {
        if (observedDataStdDev.size() == 0) {
            this->observedDataStdDev.clear();
        } else {
            throw AmiException("Input observedDataStdDev did not match dimensions nt (%i) x nytrue (%i), was %i", nt(), nytrue, observedDataStdDev.size());
        }
    }
}
    
void ExpData::setObservedDataStdDev(const realtype StdDev) {
    std::fill(observedDataStdDev.begin() ,observedDataStdDev.end(), StdDev);
}
    
void ExpData::setObservedDataStdDev(const std::vector<realtype> &observedDataStdDev, int iy) {
    if (observedDataStdDev.size() == (unsigned) nt()) {
        for (int it = 0; it < nt(); ++it)
            this->observedDataStdDev.at(iy + it*nytrue) = observedDataStdDev.at(iy);
    } else {
        throw AmiException("Input observedDataStdDev did not match dimensions nt (%i), was %i", nt(), observedDataStdDev.size());
    }
}
    
void ExpData::setObservedDataStdDev(const realtype StdDev, int iy) {
    for (int it = 0; it < nt(); ++it)
        observedDataStdDev.at(iy + it*nytrue) = StdDev;
}
    
bool ExpData::isSetObservedDataStdDev(int it, int iy) const {
    if (!observedDataStdDev.empty()) {
        return !isNaN(observedDataStdDev.at(it * nytrue + iy));
    } else {
        return false;
    }
}
    
std::vector<realtype> ExpData::getObservedDataStdDev() const {
    return observedDataStdDev;
}

const realtype *ExpData::getObservedDataStdDevPtr(int it) const {
    if (!observedDataStdDev.empty()) {
        return &observedDataStdDev.at(it*nytrue);
    } else {
        return nullptr;
    }
}
    
void ExpData::setObservedEvents(const std::vector<realtype> &observedEvents) {
    if (observedEvents.size() == (unsigned) nmaxevent*nztrue) {
        this->observedEvents = observedEvents;
    } else {
        if (observedEvents.size() == 0) {
            this->observedEvents.clear();
        } else {
            throw AmiException("Input observedEvents did not match dimensions nmaxevent (%i) x nztrue (%i), was %i", nmaxevent, nztrue, observedEvents.size());
        }
    }
}
    
void ExpData::setObservedEvents(const std::vector<realtype> &observedEvents, int iz) {
    if (observedEvents.size() == (unsigned) nmaxevent) {
        for (int ie = 0; ie < nmaxevent; ++ie)
            this->observedEvents.at(iz + ie*nztrue) = observedEvents.at(iz);
    } else {
        throw AmiException("Input observedEvents did not match dimensions nmaxevent (%i), was %i", nmaxevent, observedEvents.size());
    }
}
    
bool ExpData::isSetObservedEvents(int ie, int iz) const {
    if (!observedEvents.size()) {
        return !isNaN(observedEvents.at(ie * nztrue + iz));
    } else {
        return false;
    }
}

std::vector<realtype> ExpData::getObservedEvents() const {
    return observedEvents;
}

const realtype *ExpData::getObservedEventsPtr(int ie) const {
    if (!observedEvents.empty())
        return &observedEvents.at(ie*nztrue);
    else
        return nullptr;
}
    
void ExpData::setObservedEventsStdDev(const std::vector<realtype> &observedEventsStdDev) {
    if (observedEventsStdDev.size() == (unsigned) nmaxevent*nztrue) {
        this->observedEventsStdDev = observedEventsStdDev;
    } else {
        if (observedEventsStdDev.size() == 0) {
            this->observedEventsStdDev.clear();
        } else {
            throw AmiException("Input observedEventsStdDev did not match dimensions nmaxevent (%i) x nztrue (%i), was %i", nmaxevent, nztrue, observedEventsStdDev.size());
        }
    }
}
    
void ExpData::setObservedEventsStdDev(const realtype StdDev) {
    std::fill(observedEventsStdDev.begin() ,observedEventsStdDev.end(), StdDev);
}

void ExpData::setObservedEventsStdDev(const std::vector<realtype> &observedEventsStdDev, int iz) {
    if (observedEventsStdDev.size() == (unsigned) nmaxevent) {
        for (int ie = 0; ie < nmaxevent; ++ie)
            this->observedEventsStdDev.at(iz + ie*nztrue) = observedEventsStdDev.at(iz);
    } else {
        throw AmiException("Input observedEventsStdDev did not match dimensions nmaxevent (%i), was %i", nmaxevent, observedEventsStdDev.size());
    }
}

void ExpData::setObservedEventsStdDev(const realtype StdDev, int iz) {
    for (int ie = 0; ie < nmaxevent; ++ie)
        observedEventsStdDev.at(iz + ie*nztrue) = StdDev;
}
    
bool ExpData::isSetObservedEventsStdDev(int ie, int iz) const {
    if (!observedEventsStdDev.empty()) {
        return !isNaN(observedEventsStdDev.at(ie * nztrue + iz));
    } else {
        return false;
    }
}

std::vector<realtype> ExpData::getObservedEventsStdDev() const {
    return observedEventsStdDev;
}

const realtype *ExpData::getObservedEventsStdDevPtr(int ie) const {
    if (!observedEventsStdDev.empty())
        return &observedEventsStdDev.at(ie*nztrue);
    else
        return nullptr;
}

} // namespace amici
