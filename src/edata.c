#include "include/edata.h"
#include<include/edata_accessors.h>

void freeExpData(ExpData edata) {
    if(edata) {
        free(my);
        free(ysigma);
        if(mz) free(mz);
        if(zsigma) free(zsigma);
        free(edata);
    }
}
