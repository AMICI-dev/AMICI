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


static void wrap_ErrHandlerFn(int error_code, const char *module, const char *function, char *msg, void *eh_data) {
    char buffer [250];
    sprintf(buffer,"AMICI ERROR: in module %s in function %s : %s ",module,function,msg);
    mexWarnMsgTxt(buffer);
};

static void *AMICreate(int lmm, int iter) {
    return CVodeCreate(lmm,iter);
};

static int AMISStolerances(void *mem, double rtol,double atol) {
    return CVodeSStolerances(mem,rtol,atol);
};

static int AMISensEEtolerances(void *mem) {
    return CVodeSensEEtolerances(mem);
};

static int AMISetSensErrCon(void *mem,bool error_corr) {
    return CVodeSetSensErrCon(mem,error_corr);
};

static int AMISetQuadErrConB(void *mem,int which, bool flag) {
    return CVodeSetQuadErrConB(mem,which,flag);
};

static int AMIGetRootInfo(void *mem,int *rootsfound) {
    return CVodeGetRootInfo(mem,rootsfound);
};

static int AMISetErrHandlerFn(void *mem) {
    return CVodeSetErrHandlerFn(mem,wrap_ErrHandlerFn,NULL);
};

static int AMISetUserData(void *mem, void *user_data) {
    return CVodeSetUserData(mem,user_data);
};

static int AMISetUserDataB(void *mem, int which, void *user_data) {
    return CVodeSetUserDataB(mem, which, user_data);
};

static int AMISetMaxNumSteps(void *mem, long int mxsteps) {
    return CVodeSetMaxNumSteps(mem,mxsteps);
};

static int AMISetStabLimDet(void *mem, int stldet) {
    return CVodeSetStabLimDet(mem,stldet);
};

static int AMISetStabLimDetB(void *mem, int which, int stldet) {
    return CVodeSetStabLimDetB(mem,which,stldet);
};

static int AMISetId(void *mem, N_Vector id) {
    return(0);
};

static int AMISetSuppressAlg(void *mem, bool flag) {
    return(0);
};

static int AMIReInit(void *mem, realtype t0, N_Vector yy0, N_Vector yp0) {
    return CVodeReInit( mem, t0, yy0);
};

static int AMISensReInit(void *mem, int ism, N_Vector *yS0, N_Vector *ypS0) {
    return CVodeSensReInit( mem, ism, yS0);
};

static int AMISetSensParams(void *mem, realtype *p, realtype *pbar, int *plist) {
    return CVodeSetSensParams(mem, p, pbar, plist);
};

static int AMIGetDky(void *mem, realtype t, int k, N_Vector dky) {
    return CVodeGetDky(mem, t, k, dky);
};

static int AMIGetSens(void *mem, realtype *tret, N_Vector *yySout) {
    return CVodeGetSens( mem, tret, yySout);
};

static int AMIRootInit(void *mem, int nrtfn, void *ptr) {
    return CVodeRootInit( mem, nrtfn, ptr);
};

static void AMIFree(void **mem) {
    CVodeFree(mem);
};

static int AMIAdjInit(void *mem, long int steps, int interp) {
    return CVodeAdjInit(mem, steps, interp);
};

static int AMICreateB(void *mem, int lmm, int iter, int *which) {
    return CVodeCreateB(mem, lmm, iter, which);
};

static int AMIReInitB(void *mem, int which, realtype tB0, N_Vector yyB0, N_Vector ypB0) {
    return CVodeReInitB(mem, which, tB0, yyB0);
};

static int AMISStolerancesB(void *mem, int which, realtype relTolB, realtype absTolB) {
    return CVodeSStolerancesB(mem, which, relTolB, absTolB);
};

static int AMIQuadReInitB(void *mem, int which, N_Vector yQB0) {
    return CVodeQuadReInitB(mem, which, yQB0);
};

static int AMIQuadSStolerancesB(void *mem, int which, realtype reltolQB, realtype abstolQB) {
    return CVodeQuadSStolerancesB(mem, which, reltolQB, abstolQB);
};

static int AMISolve(void *mem, realtype tout, N_Vector yret, N_Vector ypret, realtype *tret, int itask) {
    return CVode(mem, tout, yret, tret, itask);
};

static int AMISolveF(void *mem, realtype tout, N_Vector yret, N_Vector ypret, realtype *tret, int itask, int *ncheckPtr) {
    return CVodeF(mem, tout, yret, tret, itask, ncheckPtr);
};

static int AMISolveB(void *mem, realtype tBout, int itaskB) {
    return CVodeB(mem, tBout, itaskB);
};

static int AMISetMaxNumStepsB(void *mem, int which, long int mxstepsB) {
    return CVodeSetMaxNumStepsB(mem, which, mxstepsB);
};

static int AMIGetB(void *mem, int which, realtype *tret, N_Vector yy, N_Vector yp) {
    return CVodeGetB(mem, which, tret, yy);
};

static int AMIGetQuadB(void *mem, int which, realtype *tret, N_Vector qB) {
    return CVodeGetQuadB(mem, which, tret, qB);
};

static int AMIDense(void *mem, int nx) {
    return CVDense(mem, nx);
};

static int AMIDenseB(void *mem, int which, int nx) {
    return CVDenseB(mem, which, nx);
};

static int AMIBand(void *mem, int nx, int ubw, int lbw) {
    return CVBand(mem, nx, ubw, lbw);
};

static int AMIBandB(void *mem, int which, int nx, int ubw, int lbw) {
    return CVBandB(mem, which, nx, ubw, lbw);
};

static int AMIDiag(void *mem) {
    return CVDiag(mem);
};

static int AMIDiagB(void *mem, int which) {
    return CVDiagB(mem,which);
};

static int AMISpgmr(void *mem, int prectype, int maxl) {
    return CVSpgmr(mem, prectype, maxl);
};

static int AMISpgmrB(void *mem, int which, int prectype, int maxl) {
    return CVSpgmrB(mem, which, prectype, maxl);
};

static int AMISpbcg(void *mem, int prectype, int maxl) {
    return CVSpbcg(mem, prectype, maxl);
};

static int AMISpbcgB(void *mem, int which, int prectype, int maxl) {
    return CVSpbcgB(mem, which, prectype, maxl);
};

static int AMISptfqmr(void *mem, int prectype, int maxl) {
    return CVSptfqmr(mem, prectype, maxl);
};

static int AMISptfqmrB(void *mem, int which, int prectype, int maxl) {
    return CVSptfqmrB(mem, which, prectype, maxl);
};

static int AMIKLU(void *mem, int nx, int nnz) {
    return CVKLU(mem, nx, nnz);
};

static int AMIKLUSetOrdering(void *mem, int ordering) {
    return CVKLUSetOrdering(mem, ordering);
};

static int AMIKLUSetOrderingB(void *mem, int which, int ordering) {
    return CVKLUSetOrderingB(mem, which, ordering);
};

static int AMIKLUB(void *mem, int which, int nx, int nnz) {
    return CVKLUB(mem, which, nx, nnz);
};

static int AMIGetNumSteps(void *mem, long int *numsteps) {
    return CVodeGetNumSteps(mem,numsteps);
}

static int AMIGetNumRhsEvals(void *mem, long int *numrhsevals) {
    return CVodeGetNumRhsEvals(mem,numrhsevals);
}

static int AMIGetLastOrder(void *mem,int *order) {
    return CVodeGetLastOrder(mem,order);
}

static void *AMIGetAdjBmem(void *mem, int which) {
    return CVodeGetAdjCVodeBmem(mem,which);
}

static int AMICalcIC(void *mem, realtype tout1) {
    return(0);
}

static int AMICalcICB(void *mem, int which, realtype tout1, N_Vector xB, N_Vector dxB) {
    return(0);
}

#endif /* CVodewrap_h */