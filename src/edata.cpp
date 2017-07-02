#include "include/edata.h"

#include <include/udata.h>
#include <include/amici.h>



ExpData::ExpData(const UserData *udata) {
    my = sigmay = mz = sigmaz = NULL;
    if(udata){
        my = new double[udata->nt*udata->nytrue]();
        sigmay = new double[udata->nt*udata->nytrue]();
        mz = new double[udata->nmaxevent*udata->nztrue]();
        mrz = new double[udata->nmaxevent*udata->nztrue]();
        sigmaz = new double[udata->nmaxevent*udata->nztrue]();
    }
}

ExpData::~ExpData() {
    if(my) delete[] my;
    if(sigmay) delete[] sigmay;
    if(mz) delete[] mz;
    if(mrz) delete[] mrz;
    if(sigmaz) delete[] sigmaz;
}
