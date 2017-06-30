#include "include/edata.h"

#include <include/udata.h>
#include <include/amici.h>



ExpData::ExpData(const UserData *udata) {
    my = new double[udata->nt*udata->ny]();
    sigmay = new double[udata->nt*udata->ny]();
    mz = new double[udata->nmaxevent*udata->nz]();
    sigmaz = new double[udata->nmaxevent*udata->nz]();
}

ExpData::~ExpData() {
    if(my) delete[] my;
    if(sigmay) delete[] sigmay;
    if(mz) delete[] mz;
    if(sigmaz) delete[] sigmaz;
}
