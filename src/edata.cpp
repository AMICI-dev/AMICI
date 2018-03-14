#include "amici/edata.h"

#include "amici/defines.h"
#include "amici/model.h"

#include <cstring>

namespace amici {

ExpData::ExpData() : nytrue(0), nztrue(0), nt(0), nmaxevent(0) {}

ExpData::ExpData(Model const& model)
    : nytrue(model.nytrue),
      nztrue(model.nztrue),
      nt(model.nt()),
      nmaxevent(model.nMaxEvent()) {
    /**
     * constructor that initializes with Model
     *
     * @param model pointer to model specification object @type Model
     */
    my.resize(model.nt() * model.nytrue);
    sigmay.resize(model.nt() * model.nytrue);
    mz.resize(model.nMaxEvent() * model.nztrue);
    sigmaz.resize(model.nMaxEvent() * model.nztrue);
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
