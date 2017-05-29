#include "include/amici_interface_cpp.h"

#include <include/edata_accessors.h>
#include <include/rdata_accessors.h>
#include <include/tdata_accessors.h>

#include <cstring>

/**
 * @ brief initialise matrix and attach to the field
 * @ param FIELD name of the field to which the matrix will be attached
 * @ param D1 number of rows in the matrix
 * @ param D2 number of columns in the matrix
 */
#define initField2(FIELD,D1,D2) \
FIELD ## data = new double[(D1) * (D2)]();

/**
 * @ brief initialise 3D tensor and attach to the field
 * @ param FIELD name of the field to which the tensor will be attached
 * @ param D1 number of rows in the tensor
 * @ param D2 number of columns in the tensor
 * @ param D3 number of elements in the third dimension of the tensor
 */
#define initField3(FIELD,D1,D2,D3) \
FIELD ## data = new double[(D1) * (D2) * (D3)]();

/**
 * @ brief initialise 4D tensor and attach to the field
 * @ param FIELD name of the field to which the tensor will be attached
 * @ param D1 number of rows in the tensor
 * @ param D2 number of columns in the tensor
 * @ param D3 number of elements in the third dimension of the tensor
 * @ param D4 number of elements in the fourth dimension of the tensor
 */
#define initField4(FIELD,D1,D2,D3,D4) \
FIELD ## data = new double[(D1) * (D2) * (D3) * (D4)]();

ReturnData *initReturnData(const UserData *udata, int *pstatus) {
    ReturnData *rdata; /* returned rdata struct */

    /* Return rdata structure */
    rdata = new ReturnData();
    if (rdata == NULL)
        return(NULL);

    memset(rdata, 0, sizeof(*rdata));

    tsdata = new double[udata->nt]();

    #include "include/amici_init_return_data_fields.h"

    return(rdata);
}


ReturnData *getSimulationResults(UserData *udata, const ExpData *edata, int *pstatus) {
    double *originalParams = NULL;

    if(udata->pscale != AMI_SCALING_NONE) {
        originalParams = (double *) malloc(sizeof(double) * udata->np);
        memcpy(originalParams, udata->p, sizeof(double) * udata->np);
    }

    ReturnData *rdata = initReturnData(udata, pstatus);

    runAmiciSimulation(udata, edata, rdata, pstatus);

    if(originalParams) {
        memcpy(udata->p, originalParams, sizeof(double) * udata->np);
        free(originalParams);
    }

    return rdata;
}
