#include "include/edata.h"
#include<include/edata_accessors.h>

void freeExpData(ExpData *edata) {
    if(edata) {
        if(my) delete my;
        if(ysigma) delete ysigma;
        if(mz) delete mz;
        if(zsigma) delete zsigma;
        delete edata;
    }
}
