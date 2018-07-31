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
    : ts(std::move(ts)), my(std::move(my)), sigmay(std::move(sigmay)), mz(std::move(mz)), sigmaz(std::move(sigmaz)),
      nytrue(nytrue), nztrue(nztrue), nmaxevent(nmaxevent)
{
}

ExpData::ExpData(Model const& model)
    : ExpData(model.nytrue, model.nztrue, model.nMaxEvent(), model.getTimepoints())
{
    fixedParameters = model.getFixedParameters();
}

ExpData::ExpData(const ExpData &other)
    : ExpData(other.nytrue, other.nztrue, other.nmaxevent,
              other.ts, other.my, other.sigmay, other.mz, other.sigmaz)
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
    if (sigma_y.size() != nytrue)
        throw AmiException("Dimension of sigma_y must be %d, was %d", nytrue, sigma_y.size());
        
    if (sigma_z.size() != nztrue)
        throw AmiException("Dimension of sigma_z must be %d, was %d", nztrue, sigma_z.size());
    
    std::random_device rd{};
    std::mt19937 gen{rd()};
    
    for (int iy = 0; iy < nytrue; ++iy) {
        if (sigma_y[iy] < 0.0)
            throw AmiException("All sigma_y must positive! sigma at index %d was %f", iy, sigma_y[iy]);
        std::normal_distribution<> e{0, sigma_y[iy]};
        for (int it = 0; it < nt(); ++it) {
            my.at(iy + rdata.nytrue * it) = rdata.y.at(iy + rdata.ny * it) + e(gen);
            sigmay.at(iy + rdata.ny * it) = sigma_y[iy];
        }
    }
    
    for (int iz = 0; iz < nztrue; ++iz) {
        if (sigma_z[iz] < 0.0)
            throw AmiException("All sigma_z must positive! sigma at index %d was %f", iz, sigma_z[iz]);
        std::normal_distribution<> e{0, sigma_z[iz]};
        for (int ie = 0; ie < nmaxevent; ++ie) {
            mz.at(iz + rdata.nztrue * ie) = rdata.z.at(iz + rdata.nz * ie) + e(gen);
            sigmay.at(iz + rdata.nztrue * ie) = sigma_z[iz];
        }
    }
}

void ExpData::setTimepoints(const std::vector<realtype> &ts) {
    this->ts = ts;
    my.clear();
    my.resize(nt()*nytrue, getNaN());
    sigmay.clear();
    sigmay.resize(nt()*nytrue, getNaN());
}
    
std::vector<realtype> ExpData::getTimepoints() const {
    return ts;
}
    
realtype ExpData::getTimepoint(int it) const {
    return ts.at(it);
}
    
void ExpData::setObservedData(const std::vector<realtype> &observedData) {
    if (observedData.size() == nt()*nytrue) {
        my = observedData;
    } else {
        if (observedData.size() == 0) {
            my.clear();
        } else {
            throw AmiException("Input observedData did not match dimensions nt (%i) x nytrue (%i), was %i", nt(), nytrue, observedData.size());
        }
    }
}
    
void ExpData::setObservedData(const std::vector<realtype> &observedData, int iy) {
    if (observedData.size() == nt()) {
        for (int it = 0; it < nt(); ++it)
            my.at(iy + it*nytrue) = observedData.at(iy);
    } else {
        throw AmiException("Input observedData did not match dimensions nt (%i), was %i", nt(), observedData.size());
    }
}
    
std::vector<realtype> ExpData::getObservedData() const {
    return my;
}
    
const realtype *ExpData::getObservedData(int it) const {
    return &my.at(it*nytrue);
}

void ExpData::setObservedDataStdDev(const std::vector<realtype> &observedDataStdDev) {
    if (observedDataStdDev.size() == nt()*nytrue) {
        sigmay = observedDataStdDev;
    } else {
        if (observedDataStdDev.size() == 0) {
            sigmay.clear();
        } else {
            throw AmiException("Input observedDataStdDev did not match dimensions nt (%i) x nytrue (%i), was %i", nt(), nytrue, observedDataStdDev.size());
        }
    }
}
    
void ExpData::setObservedDataStdDev(const realtype StdDev) {
    std::fill(sigmay.begin() ,sigmay.end(), StdDev);
}
    
void ExpData::setObservedDataStdDev(const std::vector<realtype> &observedDataStdDev, int iy) {
    if (observedDataStdDev.size() == nt()) {
        for (int it = 0; it < nt(); ++it)
            sigmay.at(iy + it*nytrue) = observedDataStdDev.at(iy);
    } else {
        throw AmiException("Input observedDataStdDev did not match dimensions nt (%i), was %i", nt(), observedDataStdDev.size());
    }
}
    
void ExpData::setObservedDataStdDev(const realtype StdDev, int iy) {
    for (int it = 0; it < nt(); ++it)
        sigmay.at(iy + it*nytrue) = StdDev;
}
    
std::vector<realtype> ExpData::getObservedDataStdDev() const {
    return sigmay;
}

const realtype *ExpData::getObservedDataStdDev(int it) const {
    return &sigmay.at(it*nytrue);
}
    
void ExpData::setObservedEvents(const std::vector<realtype> &observedEvents) {
    if (observedEvents.size() == nmaxevent*nztrue) {
        mz = observedEvents;
    } else {
        if (observedEvents.size() == 0) {
            mz.clear();
        } else {
            throw AmiException("Input observedEvents did not match dimensions nmaxevent (%i) x nztrue (%i), was %i", nmaxevent, nztrue, observedEvents.size());
        }
    }
}
    
void ExpData::setObservedEvents(const std::vector<realtype> &observedEvents, int iz) {
    if (observedEvents.size() == nmaxevent) {
        for (int ie = 0; ie < nmaxevent; ++ie)
            mz.at(iz + ie*nztrue) = observedEvents.at(iz);
    } else {
        throw AmiException("Input observedEvents did not match dimensions nmaxevent (%i), was %i", nmaxevent, observedEvents.size());
    }
}

std::vector<realtype> ExpData::getObservedEvents() const {
    return mz;
}

const realtype *ExpData::getObservedEvents(int ie) const {
    return &mz.at(ie*nztrue);
}
    
void ExpData::setObservedEventsStdDev(const std::vector<realtype> &observedEventsStdDev) {
    if (observedEventsStdDev.size() == nmaxevent*nztrue) {
        sigmaz = observedEventsStdDev;
    } else {
        if (observedEventsStdDev.size() == 0) {
            sigmaz.clear();
        } else {
            throw AmiException("Input observedEventsStdDev did not match dimensions nmaxevent (%i) x nztrue (%i), was %i", nmaxevent, nztrue, observedEventsStdDev.size());
        }
    }
}
    
void ExpData::setObservedEventsStdDev(const realtype StdDev) {
    std::fill(sigmaz.begin() ,sigmaz.end(), StdDev);
}

void ExpData::setObservedEventsStdDev(const std::vector<realtype> &observedEventsStdDev, int iz) {
    if (observedEventsStdDev.size() == nmaxevent) {
        for (int ie = 0; ie < nmaxevent; ++ie)
            sigmaz.at(iz + ie*nztrue) = observedEventsStdDev.at(iz);
    } else {
        throw AmiException("Input observedEventsStdDev did not match dimensions nmaxevent (%i), was %i", nmaxevent, observedEventsStdDev.size());
    }
}

void ExpData::setObservedEventsStdDev(const realtype StdDev, int iz) {
    for (int ie = 0; ie < nmaxevent; ++ie)
        sigmaz.at(iz + ie*nztrue) = StdDev;
}

std::vector<realtype> ExpData::getObservedEventsStdDev() const {
    return sigmaz;
}

const realtype *ExpData::getObservedEventsStdDev(int ie) const {
    return &sigmaz.at(ie*nztrue);
}

} // namespace amici
