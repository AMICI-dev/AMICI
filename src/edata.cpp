#include "include/edata.h"

#include "include/amici_defines.h"
#include "include/amici_model.h"
#include <include/udata.h>

#include <cstring>

namespace amici {

ExpData::ExpData() : nytrue(0), nztrue(0), nt(0), nmaxevent(0) {}

ExpData::ExpData(const UserData *udata, Model *model)
    : nytrue(model->nytrue), nztrue(model->nztrue), nt(udata->nt()),
      nmaxevent(udata->nme()) {
    /**
     * constructor that initializes with UserData and model
     *
     * @param[out] udata pointer to the return data struct @type ReturnData
     * @param[in] model pointer to model specification object @type Model
     */
    if (udata) {
        my.resize(udata->nt() * model->nytrue);
        sigmay.resize(udata->nt() * model->nytrue);
        mz.resize(udata->nme() * model->nztrue);
        sigmaz.resize(udata->nme() * model->nztrue);
    }
}

void ExpData::setObservedData(const double *observedData) {
    /**
     * set function that copies data from input to ExpData::my
     *
     * @param[in] observedData observed data
     */
    for (int imy = 0; imy < nytrue * nt; ++imy) {
        my.at(imy) = static_cast<const realtype>(observedData[imy]);
    }
}

void ExpData::setObservedDataStdDev(const double *observedDataStdDev) {
    /**
     * set function that copies data from input to ExpData::sigmay
     *
     * @param[in] observedDataStdDev standard deviation of observed data
     */
    for (int imy = 0; imy < nytrue * nt; ++imy) {
        sigmay.at(imy) = static_cast<const realtype>(observedDataStdDev[imy]);
    }
}

void ExpData::setObservedEvents(const double *observedEvents) {
    /**
     * set function that copies data from input to ExpData::mz
     *
     * @param[in] observedEvents observed event data
     */
    for (int imz = 0; imz < nztrue * nmaxevent; ++imz) {
        mz.at(imz) = static_cast<const realtype>(observedEvents[imz]);
    }
}

void ExpData::setObservedEventsStdDev(const double *observedEventsStdDev) {
    /**
     * set function that copies data from input to ExpData::sigmaz
     *
     * @param[in] observedEventsStdDev standard deviation of observed event data
     */
    for (int imz = 0; imz < nztrue * nmaxevent; ++imz) {
        sigmaz.at(imz) = static_cast<const realtype>(observedEventsStdDev[imz]);
    }
}

ExpData::~ExpData() {
}

} // namespace amici
