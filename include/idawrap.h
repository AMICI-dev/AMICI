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
#include <include/amici_solver.h>
#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

class IDASolver : public Solver
{
public:
    IDASolver() {}

    void *AMICreate(int lmm, int iter) {
        return IDACreate();
    }

    int AMISStolerances(void *mem, double rtol,double atol) {
        return IDASStolerances(mem,rtol,atol);
    }

    int AMISensEEtolerances(void *mem) {
        return IDASensEEtolerances(mem);
    }

    int AMISetSensErrCon(void *mem,bool error_corr) {
        return IDASetSensErrCon(mem,error_corr);
    }

    int AMISetQuadErrConB(void *mem,int which, bool flag) {
        return IDASetQuadErrConB(mem,which,flag);
    }

    int AMIGetRootInfo(void *mem,int *rootsfound) {
        return IDAGetRootInfo(mem,rootsfound);
    }

    int AMISetErrHandlerFn(void *mem) {
        return IDASetErrHandlerFn(mem,wrap_ErrHandlerFn,NULL);
    }

    int AMISetUserData(void *mem, void *user_data) {
        return IDASetUserData(mem,user_data);
    }

    int AMISetUserDataB(void *mem, int which, void *user_data) {
        return IDASetUserDataB(mem, which, user_data);
    }

    int AMISetMaxNumSteps(void *mem, long int mxsteps) {
        return IDASetMaxNumSteps(mem,mxsteps);
    }

    int AMISetStabLimDet(void *mem, int stldet) {
        return(0);
    }

    int AMISetStabLimDetB(void *mem, int which, int stldet) {
        return(0);
    }

    int AMISetId(void *mem, N_Vector id) {
        return IDASetId(mem, id);
    }

    int AMISetSuppressAlg(void *mem, bool flag) {
        return IDASetSuppressAlg(mem, flag);
    }

    int AMIReInit(void *mem, realtype t0, N_Vector yy0, N_Vector yp0) {
        return IDAReInit( mem, t0, yy0, yp0);
    }

    int AMISensReInit(void *mem, int ism, N_Vector *yS0, N_Vector *ypS0) {
        return IDASensReInit( mem, ism, yS0, ypS0);
    }

    int AMISetSensParams(void *mem, realtype *p, realtype *pbar, int *plist) {
        return IDASetSensParams(mem, p, pbar, plist);
    }

    int AMIGetDky(void *mem, realtype t, int k, N_Vector dky) {
        return IDAGetDky(mem, t, k, dky);
    }

    int AMIGetSens(void *mem, realtype *tret, N_Vector *yySout) {
        return IDAGetSens( mem, tret, yySout);
    }

    int AMIRootInit(void *mem, int nrtfn, IDARootFn ptr) {
        return IDARootInit( mem, nrtfn, ptr);
    }

    void AMIFree(void **mem) {
        IDAFree(mem);
    }

    int AMIAdjInit(void *mem, long int steps, int interp) {
        return IDAAdjInit(mem, steps, interp);
    }

    int AMICreateB(void *mem, int lmm, int iter, int *which) {
        return IDACreateB(mem, which);
    }

    int AMIReInitB(void *mem, int which, realtype tB0, N_Vector yyB0, N_Vector ypB0) {
        return IDAReInitB(mem, which, tB0, yyB0, ypB0);
    }

    int AMISStolerancesB(void *mem, int which, realtype relTolB, realtype absTolB) {
        return IDASStolerancesB(mem, which, relTolB, absTolB);
    }

    int AMIQuadReInitB(void *mem, int which, N_Vector yQB0) {
        return IDAQuadReInitB(mem, which, yQB0);
    }

    int AMIQuadSStolerancesB(void *mem, int which, realtype reltolQB, realtype abstolQB) {
        return IDAQuadSStolerancesB(mem, which, reltolQB, abstolQB);
    }

    int AMISolve(void *mem, realtype tout, N_Vector yret, N_Vector ypret, realtype *tret, int itask) {
        return IDASolve(mem, tout, tret, yret, ypret, itask);
    }

    int AMISolveF(void *mem, realtype tout, N_Vector yret, N_Vector ypret, realtype *tret, int itask, int *ncheckPtr) {
        return IDASolveF(mem, tout, tret, yret, ypret, itask, ncheckPtr);
    }

    int AMISolveB(void *mem, realtype tBout, int itaskB) {
        return IDASolveB(mem, tBout, itaskB);
    }

    int AMISetMaxNumStepsB(void *mem, int which, long int mxstepsB) {
        return IDASetMaxNumStepsB(mem, which, mxstepsB);
    }

    int AMIGetB(void *mem, int which, realtype *tret, N_Vector yy, N_Vector yp) {
        return IDAGetB(mem, which, tret, yy, yp);
    }

    int AMIGetQuadB(void *mem, int which, realtype *tret, N_Vector qB) {
        return IDAGetQuadB(mem, which, tret, qB);
    }

    int AMIDense(void *mem, int nx) {
        return IDADense(mem, nx);
    }

    int AMIDenseB(void *mem, int which, int nx) {
        return IDADenseB(mem, which, nx);
    }

    int AMIBand(void *mem, int nx, int ubw, int lbw) {
        return IDABand(mem, nx, ubw, lbw);
    }

    int AMIBandB(void *mem, int which, int nx, int ubw, int lbw) {
        return IDABandB(mem, which, nx, ubw, lbw);
    }

    int AMIDiag(void *mem) {
        return(-99);
    }

    int AMIDiagB(void *mem, int which) {
        return(-99);
    }

    int AMISpgmr(void *mem, int prectype, int maxl) {
        return IDASpgmr(mem, maxl);
    }

    int AMISpgmrB(void *mem, int which, int prectype, int maxl) {
        return IDASpgmrB(mem, which, maxl);
    }

    int AMISpbcg(void *mem, int prectype, int maxl) {
        return IDASpbcg(mem, maxl);
    }

    int AMISpbcgB(void *mem, int which, int prectype, int maxl) {
        return IDASpbcgB(mem, which, maxl);
    }

    int AMISptfqmr(void *mem, int prectype, int maxl) {
        return IDASptfqmr(mem, maxl);
    }

    int AMISptfqmrB(void *mem, int which, int prectype, int maxl) {
        return IDASptfqmrB(mem, which, maxl);
    }

    int AMIKLU(void *mem, int nx, int nnz, int sparsetype) {
        return IDAKLU(mem, nx, nnz, sparsetype);
    }

    int AMIKLUSetOrdering(void *mem, int ordering) {
        return IDAKLUSetOrdering(mem, ordering);
    }

    int AMIKLUSetOrderingB(void *mem, int which, int ordering) {
        return IDAKLUSetOrderingB(mem, which, ordering);
    }

    int AMIKLUB(void *mem, int which, int nx, int nnz, int sparsetype) {
        return IDAKLUB(mem, which, nx, nnz, sparsetype);
    }

    int AMIGetNumSteps(void *mem, long int *numsteps) {
        return IDAGetNumSteps(mem,numsteps);
    }

    int AMIGetNumRhsEvals(void *mem, long int *numrhsevals) {
        return IDAGetNumResEvals(mem,numrhsevals);
    }

    int AMIGetNumErrTestFails(void *mem, long int *numerrtestfails) {
        return IDAGetNumErrTestFails(mem,numerrtestfails);
    }

    int AMIGetNumNonlinSolvConvFails(void *mem, long int *numnonlinsolvconvfails) {
        return IDAGetNumNonlinSolvConvFails(mem,numnonlinsolvconvfails);
    }

    int AMIGetLastOrder(void *mem,int *order) {
        return IDAGetLastOrder(mem,order);
    }

    void *AMIGetAdjBmem(void *mem, int which) {
        return IDAGetAdjIDABmem(mem,which);
    }

    int AMICalcIC(void *mem, realtype tout1){
        return IDACalcIC(mem,IDA_YA_YDP_INIT , tout1);
    }

    int AMICalcICB(void *mem, int which, realtype tout1, N_Vector xB, N_Vector dxB) {
        return IDACalcICB(mem, which, tout1, xB, dxB);
    }

    int AMISetStopTime(void *mem, realtype tstop) {
        return IDASetStopTime(mem, tstop);
    }

    ~IDASolver() {}
}

#endif /* idawrap_h */
