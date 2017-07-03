#ifndef idawrap_h
#define idawrap_h

#include <idas/idas.h>
/*#include <idas/idas_lapack.h>*/
#include <idas/idas_band.h>
#include <idas/idas_bbdpre.h>
#include <idas/idas_dense.h>
/*#include <idas/idas_diag.h>*/
#include <idas/idas_spbcgs.h>
#include <idas/idas_spgmr.h>
#include <idas/idas_sptfqmr.h>
#include <idas/idas_klu.h>
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

static void wrap_ErrHandlerFn(int error_code, const char *module, const char *function, char *msg, void *eh_data) {
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
    
    mexWarnMsgIdAndTxt(buffid,buffer);
};

static void *AMICreate(int lmm, int iter) {
    return IDACreate();
};

static int AMISStolerances(void *mem, double rtol,double atol) {
    return IDASStolerances(mem,rtol,atol);
};

static int AMISensEEtolerances(void *mem) {
    return IDASensEEtolerances(mem);
};

static int AMISetSensErrCon(void *mem,bool error_corr) {
    return IDASetSensErrCon(mem,error_corr);
};

static int AMISetQuadErrConB(void *mem,int which, bool flag) {
    return IDASetQuadErrConB(mem,which,flag);
};

static int AMIGetRootInfo(void *mem,int *rootsfound) {
    return IDAGetRootInfo(mem,rootsfound);
};

static int AMISetErrHandlerFn(void *mem) {
    return IDASetErrHandlerFn(mem,wrap_ErrHandlerFn,NULL);
};

static int AMISetUserData(void *mem, void *user_data) {
    return IDASetUserData(mem,user_data);
};

static int AMISetUserDataB(void *mem, int which, void *user_data) {
    return IDASetUserDataB(mem, which, user_data);
};

static int AMISetMaxNumSteps(void *mem, long int mxsteps) {
    return IDASetMaxNumSteps(mem,mxsteps);
};

static int AMISetStabLimDet(void *mem, int stldet) {
    return(0);
};

static int AMISetStabLimDetB(void *mem, int which, int stldet) {
    return(0);
};

static int AMISetId(void *mem, N_Vector id) {
    return IDASetId(mem, id);
};

static int AMISetSuppressAlg(void *mem, bool flag) {
    return IDASetSuppressAlg(mem, flag);
};

static int AMIReInit(void *mem, realtype t0, N_Vector yy0, N_Vector yp0) {
    return IDAReInit( mem, t0, yy0, yp0);
};

static int AMISensReInit(void *mem, int ism, N_Vector *yS0, N_Vector *ypS0) {
    return IDASensReInit( mem, ism, yS0, ypS0);
};

static int AMISetSensParams(void *mem, realtype *p, realtype *pbar, int *plist) {
    return IDASetSensParams(mem, p, pbar, plist);
};

static int AMIGetDky(void *mem, realtype t, int k, N_Vector dky) {
    return IDAGetDky(mem, t, k, dky);
};

static int AMIGetSens(void *mem, realtype *tret, N_Vector *yySout) {
    return IDAGetSens( mem, tret, yySout);
};

static int AMIRootInit(void *mem, int nrtfn, IDARootFn ptr) {
    return IDARootInit( mem, nrtfn, ptr);
};

static void AMIFree(void **mem) {
    IDAFree(mem);
};

static int AMIAdjInit(void *mem, long int steps, int interp) {
    return IDAAdjInit(mem, steps, interp);
};

static int AMICreateB(void *mem, int lmm, int iter, int *which) {
    return IDACreateB(mem, which);
};

static int AMIReInitB(void *mem, int which, realtype tB0, N_Vector yyB0, N_Vector ypB0) {
    return IDAReInitB(mem, which, tB0, yyB0, ypB0);
};

static int AMISStolerancesB(void *mem, int which, realtype relTolB, realtype absTolB) {
    return IDASStolerancesB(mem, which, relTolB, absTolB);
};

static int AMIQuadReInitB(void *mem, int which, N_Vector yQB0) {
    return IDAQuadReInitB(mem, which, yQB0);
};

static int AMIQuadSStolerancesB(void *mem, int which, realtype reltolQB, realtype abstolQB) {
    return IDAQuadSStolerancesB(mem, which, reltolQB, abstolQB);
};

static int AMISolve(void *mem, realtype tout, N_Vector yret, N_Vector ypret, realtype *tret, int itask) {
    return IDASolve(mem, tout, tret, yret, ypret, itask);
};

static int AMISolveF(void *mem, realtype tout, N_Vector yret, N_Vector ypret, realtype *tret, int itask, int *ncheckPtr) {
    return IDASolveF(mem, tout, tret, yret, ypret, itask, ncheckPtr);
};

static int AMISolveB(void *mem, realtype tBout, int itaskB) {
    return IDASolveB(mem, tBout, itaskB);
};

static int AMISetMaxNumStepsB(void *mem, int which, long int mxstepsB) {
    return IDASetMaxNumStepsB(mem, which, mxstepsB);
};

static int AMIGetB(void *mem, int which, realtype *tret, N_Vector yy, N_Vector yp) {
    return IDAGetB(mem, which, tret, yy, yp);
};

static int AMIGetQuadB(void *mem, int which, realtype *tret, N_Vector qB) {
    return IDAGetQuadB(mem, which, tret, qB);
};

static int AMIDense(void *mem, int nx) {
    return IDADense(mem, nx);
};

static int AMIDenseB(void *mem, int which, int nx) {
    return IDADenseB(mem, which, nx);
};

static int AMIBand(void *mem, int nx, int ubw, int lbw) {
    return IDABand(mem, nx, ubw, lbw);
};

static int AMIBandB(void *mem, int which, int nx, int ubw, int lbw) {
    return IDABandB(mem, which, nx, ubw, lbw);
};

static int AMIDiag(void *mem) {
    return(-99);
};

static int AMIDiagB(void *mem, int which) {
    return(-99);
};

static int AMISpgmr(void *mem, int prectype, int maxl) {
    return IDASpgmr(mem, maxl);
};

static int AMISpgmrB(void *mem, int which, int prectype, int maxl) {
    return IDASpgmrB(mem, which, maxl);
};

static int AMISpbcg(void *mem, int prectype, int maxl) {
    return IDASpbcg(mem, maxl);
};

static int AMISpbcgB(void *mem, int which, int prectype, int maxl) {
    return IDASpbcgB(mem, which, maxl);
};

static int AMISptfqmr(void *mem, int prectype, int maxl) {
    return IDASptfqmr(mem, maxl);
};

static int AMISptfqmrB(void *mem, int which, int prectype, int maxl) {
    return IDASptfqmrB(mem, which, maxl);
};

static int AMIKLU(void *mem, int nx, int nnz, int sparsetype) {
    return IDAKLU(mem, nx, nnz, sparsetype);
};

static int AMIKLUSetOrdering(void *mem, int ordering) {
    return IDAKLUSetOrdering(mem, ordering);
};

static int AMIKLUSetOrderingB(void *mem, int which, int ordering) {
    return IDAKLUSetOrderingB(mem, which, ordering);
};

static int AMIKLUB(void *mem, int which, int nx, int nnz, int sparsetype) {
    return IDAKLUB(mem, which, nx, nnz, sparsetype);
};

static int AMIGetNumSteps(void *mem, long int *numsteps) {
    return IDAGetNumSteps(mem,numsteps);
}

static int AMIGetNumRhsEvals(void *mem, long int *numrhsevals) {
    return IDAGetNumResEvals(mem,numrhsevals);
}

static int AMIGetNumErrTestFails(void *mem, long int *numerrtestfails) {
    return IDAGetNumErrTestFails(mem,numerrtestfails);
}

static int AMIGetNumNonlinSolvConvFails(void *mem, long int *numnonlinsolvconvfails) {
    return IDAGetNumNonlinSolvConvFails(mem,numnonlinsolvconvfails);
}

static int AMIGetLastOrder(void *mem,int *order) {
    return IDAGetLastOrder(mem,order);
}

static void *AMIGetAdjBmem(void *mem, int which) {
    return IDAGetAdjIDABmem(mem,which);
}

static int AMICalcIC(void *mem, realtype tout1){
    return IDACalcIC(mem,IDA_YA_YDP_INIT , tout1);
}

static int AMICalcICB(void *mem, int which, realtype tout1, N_Vector xB, N_Vector dxB) {
    return IDACalcICB(mem, which, tout1, xB, dxB);
}

static int AMISetStopTime(void *mem, realtype tstop) {
    return IDASetStopTime(mem, tstop);
}


#endif /* idawrap_h */
