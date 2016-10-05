#include "include/edata.h"
#include<include/edata_accessors.h>

void freeExpData(ExpData *edata) {
    if(edata) {
        delete[] my;
        delete[] ysigma;
        if(mz) delete[] mz;
        if(zsigma) delete[] zsigma;
        delete edata;
    }
}
