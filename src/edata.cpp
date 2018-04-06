#include "amici/edata.h"

#include "amici/defines.h"
#include "amici/model.h"

#include <cstring>

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
}

ExpData::ExpData(const ExpData &other)
    : ExpData(other.nytrue, other.nztrue, other.nt, other.nmaxevent,
              other.my, other.sigmay, other.mz, other.sigmaz)
{
}

void ExpData::setObservedData(const double *observedData) {
    /**
     * set function that copies data from input to ExpData::my
     *
     * @param observedData observed data
     */
    for (int imy = 0; imy < nytrue * nt; ++imy) {
        my.at(imy) = static_cast<const realtype>(observedData[imy]);
    }
}

void ExpData::setObservedDataStdDev(const double *observedDataStdDev) {
    /**
     * set function that copies data from input to ExpData::sigmay
     *
     * @param observedDataStdDev standard deviation of observed data
     */
    for (int imy = 0; imy < nytrue * nt; ++imy) {
        sigmay.at(imy) = static_cast<const realtype>(observedDataStdDev[imy]);
    }
}

void ExpData::setObservedEvents(const double *observedEvents) {
    /**
     * set function that copies data from input to ExpData::mz
     *
     * @param observedEvents observed event data
     */
    for (int imz = 0; imz < nztrue * nmaxevent; ++imz) {
        mz.at(imz) = static_cast<const realtype>(observedEvents[imz]);
    }
}

void ExpData::setObservedEventsStdDev(const double *observedEventsStdDev) {
    /**
     * set function that copies data from input to ExpData::sigmaz
     *
     * @param observedEventsStdDev standard deviation of observed event data
     */
    for (int imz = 0; imz < nztrue * nmaxevent; ++imz) {
        sigmaz.at(imz) = static_cast<const realtype>(observedEventsStdDev[imz]);
    }
}

ExpData::~ExpData() {
}

} // namespace amici
