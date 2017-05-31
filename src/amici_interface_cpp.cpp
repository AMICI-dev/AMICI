#include "include/amici_interface_cpp.h"

#include <include/edata_accessors.h>
#include <include/tdata_accessors.h>

#include <cstring>

#define initField2(FIELD,D1,D2) \
FIELD ## data = new double[(D1) * (D2)]();
#define initField3(FIELD,D1,D2,D3) \
FIELD ## data = new double[(D1) * (D2) * (D3)]();

#define initField4(FIELD,D1,D2,D3,D4) \
FIELD ## data = new double[(D1) * (D2) * (D3) * (D4)]();


ReturnData *getSimulationResults(UserData *udata, const ExpData *edata) {
    double *originalParams = NULL;

    if(udata->pscale != AMI_SCALING_NONE) {
        originalParams = (double *) malloc(sizeof(double) * udata->np);
        memcpy(originalParams, udata->p, sizeof(double) * udata->np);
    }

    ReturnData *rdata = new ReturnData(udata);

    int status;
    runAmiciSimulation(udata, edata, rdata, &status);
    *rdata->status = status;

    if(originalParams) {
        memcpy(udata->p, originalParams, sizeof(double) * udata->np);
        free(originalParams);
    }

    return rdata;
}
