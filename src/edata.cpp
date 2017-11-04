#include "include/edata.h"

#include "include/amici_defines.h"
#include "include/amici_model.h"
#include <include/udata.h>

#include <cstring>

namespace amici {

ExpData::ExpData() : nytrue(0), nztrue(0), nt(0), nmaxevent(0) {}

ExpData::ExpData(const UserData *udata, Model *model)
    : nytrue(model->nytrue), nztrue(model->nztrue), nt(udata->nt),
      nmaxevent(udata->nmaxevent) {
    /**
     * constructor that initializes with UserData and model
     *
     * @param[out] udata pointer to the return data struct @type ReturnData
     * @param[in] model pointer to model specification object @type Model
     */
    if (udata) {
        my.resize(udata->nt * model->nytrue);
        sigmay.resize(udata->nt * model->nytrue);
        mz.resize(udata->nmaxevent * model->nztrue);
        sigmaz.resize(udata->nmaxevent * model->nztrue);
    }
}

void ExpData::setObservedData(const double *observedData) {
    /**
     * set function that copies data from input to ExpData::my
     *
     * @param[in] observedData observed data @type *double
     */
    memcpy(my.data(), observedData, nytrue * nt * sizeof(double));
}

void ExpData::setObservedDataStdDev(const double *observedDataStdDev) {
    /**
     * set function that copies data from input to ExpData::sigmay
     *
     * @param[in] observedDataStdDev standard deviation of observed data @type *double
     */
    memcpy(sigmay.data(), observedDataStdDev, nytrue * nt * sizeof(double));
}

void ExpData::setObservedEvents(const double *observedEvents) {
    /**
     * set function that copies data from input to ExpData::mz
     *
     * @param[in] observedEvents observed event data @type *double
     */
    memcpy(mz.data(), observedEvents, nmaxevent * nztrue * sizeof(double));
}

void ExpData::setObservedEventsStdDev(const double *observedEventsStdDev) {
    /**
     * set function that copies data from input to ExpData::sigmaz
     *
     * @param[in] observedEventsStdDev standard deviation of observed event data @type *double
     */
    memcpy(sigmaz.data(), observedEventsStdDev, nmaxevent * nztrue * sizeof(double));
}

ExpData::~ExpData() {
}

} // namespace amici
