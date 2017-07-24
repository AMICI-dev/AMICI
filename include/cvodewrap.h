#ifndef cvodewrap_h
#define cvodewrap_h

#include <cvodes/cvodes.h>
/*#include <cvodes/cvodes_lapack.h>*/
#include <cvodes/cvodes_band.h>
#include <cvodes/cvodes_bbdpre.h>
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_diag.h>
#include <cvodes/cvodes_spbcgs.h>
#include <cvodes/cvodes_spgmr.h>
#include <cvodes/cvodes_sptfqmr.h>
#include <cvodes/cvodes_klu.h>
#include <klu.h>
#include <amd.h>
#include <colamd.h>
#include <btf.h>

#include <include/amici.h>

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

void wrap_ErrHandlerFn(int error_code, const char *module, const char *function, char *msg, void *eh_data) {
    char buffer[250];
    char buffid[250];
    sprintf(buffer,"AMICI ERROR: in module %s in function %s : %s ",module,function,msg);
    switch (error_code) {
        case 99:
            sprintf(buffid,"AMICI:mex:%s:%s:CV_WARNING",module,function);
            break;
            
        case -1:
            sprintf(buffid,"AMICI:mex:%s:%s:CV_TOO_MUCH_WORK",module,function);
            break;
            
        case -2:
            sprintf(buffid,"AMICI:mex:%s:%s:CV_TOO_MUCH_ACC",module,function);
            break;

        case -3:
            sprintf(buffid,"AMICI:mex:%s:%s:CV_ERR_FAILURE",module,function);
            break;

        case -4:
            sprintf(buffid,"AMICI:mex:%s:%s:CV_CONV_FAILURE",module,function);
            break;
            
        default:
            sprintf(buffid,"AMICI:mex:%s:%s:CV_OTHER",module,function);
            break;
    }
    
    warnMsgIdAndTxt(buffid,buffer);
};

void *AMICreate(int lmm, int iter) {
    return CVodeCreate(lmm,iter);
};

int AMISStolerances(void *mem, double rtol,double atol) {
    return CVodeSStolerances(mem,rtol,atol);
};

int AMISensEEtolerances(void *mem) {
    return CVodeSensEEtolerances(mem);
};

int AMISetSensErrCon(void *mem,bool error_corr) {
    return CVodeSetSensErrCon(mem,error_corr);
};

int AMISetQuadErrConB(void *mem,int which, bool flag) {
    return CVodeSetQuadErrConB(mem,which,flag);
};

int AMIGetRootInfo(void *mem,int *rootsfound) {
    return CVodeGetRootInfo(mem,rootsfound);
};

int AMISetErrHandlerFn(void *mem) {
    return CVodeSetErrHandlerFn(mem,wrap_ErrHandlerFn,NULL);
};

int AMISetUserData(void *mem, void *user_data) {
    return CVodeSetUserData(mem,user_data);
};

int AMISetUserDataB(void *mem, int which, void *user_data) {
    return CVodeSetUserDataB(mem, which, user_data);
};

int AMISetMaxNumSteps(void *mem, long int mxsteps) {
    return CVodeSetMaxNumSteps(mem,mxsteps);
};

int AMISetStabLimDet(void *mem, int stldet) {
    return CVodeSetStabLimDet(mem,stldet);
};

int AMISetStabLimDetB(void *mem, int which, int stldet) {
    return CVodeSetStabLimDetB(mem,which,stldet);
};

int AMISetId(void *mem, N_Vector id) {
    return(0);
};

int AMISetSuppressAlg(void *mem, bool flag) {
    return(0);
};

int AMIReInit(void *mem, realtype t0, N_Vector yy0, N_Vector yp0) {
    return CVodeReInit( mem, t0, yy0);
};

int AMISensReInit(void *mem, int ism, N_Vector *yS0, N_Vector *ypS0) {
    return CVodeSensReInit( mem, ism, yS0);
};

int AMISetSensParams(void *mem, realtype *p, realtype *pbar, int *plist) {
    return CVodeSetSensParams(mem, p, pbar, plist);
};

int AMIGetDky(void *mem, realtype t, int k, N_Vector dky) {
    return CVodeGetDky(mem, t, k, dky);
};

int AMIGetSens(void *mem, realtype *tret, N_Vector *yySout) {
    return CVodeGetSens( mem, tret, yySout);
};

int AMIRootInit(void *mem, int nrtfn, CVRootFn ptr) {
    return CVodeRootInit( mem, nrtfn, ptr);
};

void AMIFree(void **mem) {
    CVodeFree(mem);
};

int AMIAdjInit(void *mem, long int steps, int interp) {
    return CVodeAdjInit(mem, steps, interp);
};

int AMICreateB(void *mem, int lmm, int iter, int *which) {
    return CVodeCreateB(mem, lmm, iter, which);
};

int AMIReInitB(void *mem, int which, realtype tB0, N_Vector yyB0, N_Vector ypB0) {
    return CVodeReInitB(mem, which, tB0, yyB0);
};

int AMISStolerancesB(void *mem, int which, realtype relTolB, realtype absTolB) {
    return CVodeSStolerancesB(mem, which, relTolB, absTolB);
};

int AMIQuadReInitB(void *mem, int which, N_Vector yQB0) {
    return CVodeQuadReInitB(mem, which, yQB0);
};

int AMIQuadSStolerancesB(void *mem, int which, realtype reltolQB, realtype abstolQB) {
    return CVodeQuadSStolerancesB(mem, which, reltolQB, abstolQB);
};

int AMISolve(void *mem, realtype tout, N_Vector yret, N_Vector ypret, realtype *tret, int itask) {
    return CVode(mem, tout, yret, tret, itask);
};

int AMISolveF(void *mem, realtype tout, N_Vector yret, N_Vector ypret, realtype *tret, int itask, int *ncheckPtr) {
    return CVodeF(mem, tout, yret, tret, itask, ncheckPtr);
};

int AMISolveB(void *mem, realtype tBout, int itaskB) {
    return CVodeB(mem, tBout, itaskB);
};

int AMISetMaxNumStepsB(void *mem, int which, long int mxstepsB) {
    return CVodeSetMaxNumStepsB(mem, which, mxstepsB);
};

int AMIGetB(void *mem, int which, realtype *tret, N_Vector yy, N_Vector yp) {
    return CVodeGetB(mem, which, tret, yy);
};

int AMIGetQuadB(void *mem, int which, realtype *tret, N_Vector qB) {
    return CVodeGetQuadB(mem, which, tret, qB);
};

int AMIDense(void *mem, int nx) {
    return CVDense(mem, nx);
};

int AMIDenseB(void *mem, int which, int nx) {
    return CVDenseB(mem, which, nx);
};

int AMIBand(void *mem, int nx, int ubw, int lbw) {
    return CVBand(mem, nx, ubw, lbw);
};

int AMIBandB(void *mem, int which, int nx, int ubw, int lbw) {
    return CVBandB(mem, which, nx, ubw, lbw);
};

int AMIDiag(void *mem) {
    return CVDiag(mem);
};

int AMIDiagB(void *mem, int which) {
    return CVDiagB(mem,which);
};

int AMISpgmr(void *mem, int prectype, int maxl) {
    return CVSpgmr(mem, prectype, maxl);
};

int AMISpgmrB(void *mem, int which, int prectype, int maxl) {
    return CVSpgmrB(mem, which, prectype, maxl);
};

int AMISpbcg(void *mem, int prectype, int maxl) {
    return CVSpbcg(mem, prectype, maxl);
};

int AMISpbcgB(void *mem, int which, int prectype, int maxl) {
    return CVSpbcgB(mem, which, prectype, maxl);
};

int AMISptfqmr(void *mem, int prectype, int maxl) {
    return CVSptfqmr(mem, prectype, maxl);
};

int AMISptfqmrB(void *mem, int which, int prectype, int maxl) {
    return CVSptfqmrB(mem, which, prectype, maxl);
};

int AMIKLU(void *mem, int nx, int nnz, int sparsetype) {
    return CVKLU(mem, nx, nnz, sparsetype);
};

int AMIKLUSetOrdering(void *mem, int ordering) {
    return CVKLUSetOrdering(mem, ordering);
};

int AMIKLUSetOrderingB(void *mem, int which, int ordering) {
    return CVKLUSetOrderingB(mem, which, ordering);
};

int AMIKLUB(void *mem, int which, int nx, int nnz, int sparsetype) {
    return CVKLUB(mem, which, nx, nnz, sparsetype);
};

int AMIGetNumSteps(void *mem, long int *numsteps) {
    return CVodeGetNumSteps(mem,numsteps);
}

int AMIGetNumRhsEvals(void *mem, long int *numrhsevals) {
    return CVodeGetNumRhsEvals(mem,numrhsevals);
}

int AMIGetNumErrTestFails(void *mem, long int *numerrtestfails) {
    return CVodeGetNumErrTestFails(mem,numerrtestfails);
}

int AMIGetNumNonlinSolvConvFails(void *mem, long int *numnonlinsolvconvfails) {
    return CVodeGetNumNonlinSolvConvFails(mem,numnonlinsolvconvfails);
}

int AMIGetLastOrder(void *mem,int *order) {
    return CVodeGetLastOrder(mem,order);
}

void *AMIGetAdjBmem(void *mem, int which) {
    return CVodeGetAdjCVodeBmem(mem,which);
}

int AMICalcIC(void *mem, realtype tout1) {
    return(0);
}

int AMICalcICB(void *mem, int which, realtype tout1, N_Vector xB, N_Vector dxB) {
    return(0);
}

int AMISetStopTime(void *mem, realtype tstop) {
    return CVodeSetStopTime(mem, tstop);
}

#endif /* CVodewrap_h */
