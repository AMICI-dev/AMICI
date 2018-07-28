#include "amici/edata.h"
#include "amici/rdata.h"

#include "amici/defines.h"
#include "amici/model.h"

#include <cstring>
#include <random>

namespace amici {

ExpData::ExpData() : nytrue(0), nztrue(0), nt(0), nmaxevent(0) {}

ExpData::ExpData(int nytrue, int nztrue, int nt, int nmaxevent)
    : ExpData(nytrue, nztrue, nt, nmaxevent,
              std::vector<realtype>(nt * nytrue),
              std::vector<realtype>(nt * nytrue),
              std::vector<realtype>(nmaxevent * nztrue),
              std::vector<realtype>(nmaxevent * nztrue))

{
}

ExpData::ExpData(int nytrue, int nztrue, int nt, int nmaxevent,
                 const std::vector<realtype> &ts,
                 const std::vector<realtype> &my,
                 const std::vector<realtype> &sigmay,
                 const std::vector<realtype> &mz,
                 const std::vector<realtype> &sigmaz)
    : ts(ts), my(my), sigmay(sigmay), mz(mz), sigmaz(sigmaz),
      nytrue(nytrue), nztrue(nztrue), nt(nt), nmaxevent(nmaxevent)
{
}
    
ExpData::ExpData(int nytrue, int nztrue, int nt, int nmaxevent,
                 const std::vector<realtype> &my,
                 const std::vector<realtype> &sigmay,
                 const std::vector<realtype> &mz,
                 const std::vector<realtype> &sigmaz)
: my(my), sigmay(sigmay), mz(mz), sigmaz(sigmaz),
nytrue(nytrue), nztrue(nztrue), nt(nt), nmaxevent(nmaxevent)
{
}

ExpData::ExpData(Model const& model)
    : ExpData(model.nytrue, model.nztrue, model.nt(), model.nMaxEvent())
{
    ts = model.getTimepoints();
}

ExpData::ExpData(const ExpData &other)
    : ExpData(other.nytrue, other.nztrue, other.nt, other.nmaxevent,
              other.ts, other.my, other.sigmay, other.mz, other.sigmaz)
{
}
    
ExpData::ExpData(ReturnData const& rdata, realtype sigma_y, realtype sigma_z)
    : ExpData(rdata, std::vector<realtype>(rdata.nytrue, sigma_y), std::vector<realtype>(rdata.nztrue, sigma_z))
{
}
    
ExpData::ExpData(ReturnData const& rdata, std::vector<realtype> sigma_y, std::vector<realtype> sigma_z)
    : ExpData(rdata.nytrue, rdata.nztrue, rdata.nt, rdata.nmaxevent)
{
    if (sigma_y.size() != nytrue)
        throw AmiException("Dimension of sigma_y must be %d, was %d", nytrue, sigma_y.size());
        
    if (sigma_z.size() != nztrue)
        throw AmiException("Dimension of sigma_z must be %d, was %d", nztrue, sigma_z.size());
        
    ts = rdata.ts;
    
    std::random_device rd{};
    std::mt19937 gen{rd()};
    
    for (int iy = 0; iy < nytrue; ++iy) {
        if (sigma_y[iy] < 0.0)
            throw AmiException("All sigma_y must positive! sigma at index %d was %f", iy, sigma_y[iy]);
        std::normal_distribution<> e{0, sigma_y[iy]};
        for (int it = 0; it < nt; ++it) {
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

void ExpData::setTimepoints(const double *timepoints) {
    for (int it = 0; it < nt; ++it) {
        ts.at(it) = static_cast<const realtype>(timepoints[it]);
    }
}
    
void ExpData::setObservedData(const double *observedData) {
    for (int imy = 0; imy < nytrue * nt; ++imy) {
        my.at(imy) = static_cast<const realtype>(observedData[imy]);
    }
}

void ExpData::setObservedDataStdDev(const double *observedDataStdDev) {
    for (int imy = 0; imy < nytrue * nt; ++imy) {
        sigmay.at(imy) = static_cast<const realtype>(observedDataStdDev[imy]);
    }
}

void ExpData::setObservedEvents(const double *observedEvents) {
    for (int imz = 0; imz < nztrue * nmaxevent; ++imz) {
        mz.at(imz) = static_cast<const realtype>(observedEvents[imz]);
    }
}

void ExpData::setObservedEventsStdDev(const double *observedEventsStdDev) {
    for (int imz = 0; imz < nztrue * nmaxevent; ++imz) {
        sigmaz.at(imz) = static_cast<const realtype>(observedEventsStdDev[imz]);
    }
}

ExpData::~ExpData() {
}

} // namespace amici
