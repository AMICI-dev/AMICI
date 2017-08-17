#include "include/edata.h"

#include <include/udata.h>
#include "include/amici_defines.h"
#include "include/amici_model.h"

ExpData::ExpData() {
    setDefaults();
}

ExpData::ExpData(const UserData *udata, Model *model) {
    setDefaults();
    if(udata){
        my = new double[udata->nt*model->nytrue]();
        sigmay = new double[udata->nt*model->nytrue]();
        mz = new double[udata->nmaxevent*model->nztrue]();
        mrz = new double[udata->nmaxevent*model->nztrue]();
        sigmaz = new double[udata->nmaxevent*model->nztrue]();
    }
}

void ExpData::setDefaults()
{
    my = sigmay = mz = sigmaz = nullptr;
}

ExpData::~ExpData() {
    if(my) delete[] my;
    if(sigmay) delete[] sigmay;
    if(mz) delete[] mz;
    if(mrz) delete[] mrz;
    if(sigmaz) delete[] sigmaz;
}
