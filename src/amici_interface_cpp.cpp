#include "include/amici_interface_cpp.h"

#include <include/edata_accessors.h>
#include <include/udata_accessors.h>
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
FIELD ## data = new double[D1 * D2]();

/**
 * @ brief initialise 3D tensor and attach to the field
 * @ param FIELD name of the field to which the tensor will be attached
 * @ param D1 number of rows in the tensor
 * @ param D2 number of columns in the tensor
 * @ param D3 number of elements in the third dimension of the tensor
 */
#define initField3(FIELD,D1,D2,D3) \
FIELD ## data = new double[D1 * D2 * D3]();

/**
 * @ brief initialise 4D tensor and attach to the field
 * @ param FIELD name of the field to which the tensor will be attached
 * @ param D1 number of rows in the tensor
 * @ param D2 number of columns in the tensor
 * @ param D3 number of elements in the third dimension of the tensor
 * @ param D4 number of elements in the fourth dimension of the tensor
 */
#define initField4(FIELD,D1,D2,D3,D4) \
FIELD ## data = new double[D1 * D2 * D3 * D4]();


static void initUserDataFields(const UserData user_data, ReturnData *rdata, double *pstatus);

void initUserDataFields(const UserData *udata, ReturnData *rdata) {
    initField2(llh,1,1);
    initField2(chi2,1,1);

    tsdata = new double[nt]();

    initField2(numsteps,nt,1);
    initField2(numrhsevals,nt,1);
    initField2(order,nt,1);
    if(sensi >= AMI_SENSI_ORDER_FIRST){
        initField2(numstepsS,nt,1);
        initField2(numrhsevalsS,nt,1);
    }
    if((nz>0) & (ne>0)){
        initField2(z,nmaxevent,nz);
        initField2(rz,nmaxevent,nz);
        initField2(sigmaz,nmaxevent,nz);
    }
    if(nx>0) {
        initField2(x,nt,nx);
        initField2(xdot,1,nx);
        initField2(J,nx,nx);
    }
    if(ny>0) {
        initField2(y,nt,ny);
        initField2(sigmay,nt,ny);
        if (sensi_meth == AMI_SENSI_SS) {
            initField2(dydp,ny,nplist);
            initField2(dydx,ny,nx);
            initField2(dxdotdp,nx,nplist);
        }
    }
    if(sensi >= AMI_SENSI_ORDER_FIRST) {
        initField2(sllh,nplist,1);
        if (sensi_meth == AMI_SENSI_FSA) {
            initField3(sx,nt,nx,nplist);
            if(ny>0) {
                initField3(sy,nt,ny,nplist);
                initField3(ssigmay,nt,ny,nplist);
            }
            if((nz>0) & (ne>0)){
                initField3(srz,nmaxevent,nz,nplist);
                if(sensi >= AMI_SENSI_ORDER_SECOND){
                    initField4(s2rz,nmaxevent,nz,nplist,nplist);
                }
                initField3(sz,nmaxevent,nz,nplist);
                initField3(ssigmaz,nmaxevent,nz,nplist);
            }
        }
        if (sensi_meth == AMI_SENSI_ASA) {
            if(ny>0) {
                initField3(ssigmay,nt,ny,nplist);
            }
            if((nz>0) & (ne>0)){
                initField3(ssigmaz,nmaxevent,nz,nplist);
            }
        }
        if(sensi >= AMI_SENSI_ORDER_SECOND) {
            initField2(s2llh,ng-1,nplist);
        }
    }
}


ReturnData *initReturnData(const UserData *udata, int *pstatus) {
    ReturnData *rdata; /* returned rdata struct */

    /* Return rdata structure */
    rdata = new ReturnData();
    if (rdata == NULL)
        return(NULL);

    memset(rdata, 0, sizeof(*rdata));

    double dblstatus;
    initUserDataFields(udata, rdata);
    *pstatus = (int) dblstatus;

    return(rdata);
}


ReturnData *getSimulationResults(UserData *udata, const ExpData *edata, int *pstatus) {
    double *originalParams = 0;

    if(udata->am_pscale != AMI_SCALING_NONE) {
        originalParams = (double *) malloc(sizeof(double) * np);
        memcpy(originalParams, p, sizeof(double) * np);

        unscaleParameters(udata);
    }

    int iroot = 0;
    booleantype setupBdone = false;
    *pstatus = 0;
    int problem;
    ReturnData *rdata;
    TempData *tdata = new TempData();
    void *ami_mem = 0; /* pointer to cvodes memory block */
    if (tdata == NULL) goto freturn;


    if (nx>0) {
        ami_mem = setupAMI(pstatus, udata, tdata);
        if (ami_mem == NULL) goto freturn;
    }

    rdata = initReturnData(udata, pstatus);
    if (rdata == NULL) goto freturn;

    *pstatus = 0;

    problem = workForwardProblem(udata, tdata, rdata, edata, pstatus, ami_mem, &iroot);
    if(problem)
        goto freturn;


    problem = workBackwardProblem(udata, tdata, rdata, edata, pstatus, ami_mem, &iroot, &setupBdone);
    if(problem)
        goto freturn;

    applyChainRuleFactorToSimulationResults(udata, rdata, edata);

freturn:
    freeTempDataAmiMem(udata, tdata, ami_mem, setupBdone, *pstatus);

    if(originalParams) {
        memcpy(p, originalParams, sizeof(double) * np);
        free(originalParams);
    }
    return rdata;
}
