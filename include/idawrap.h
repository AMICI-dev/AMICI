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
#include <include/amici_model_functions.h>
#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

class IDASolver : public Solver
{
public:
    IDASolver() {}

    int wrap_init(N_Vector x, N_Vector dx, realtype t){
       return IDAInit(ami_mem, resultFunction, RCONST(t), x, dx);
    }

    int wrap_binit(int which, N_Vector xB, N_Vector dxB, realtype t) {
        return IDAInitB(ami_mem, which, resultFunctionB, RCONST(t), xB, dxB);
    }

    int wrap_qbinit(int which, N_Vector qBdot) {
        return IDAQuadInitB(ami_mem, which, fqBdot, qBdot);
    }

    int wrap_RootInit(UserData *udata) {
        return IDARootInit(ami_mem, udata->ne, froot);
    }

    int wrap_SensInit1(N_Vector *sx, N_Vector *sdx, UserData *udata) {
        return IDASensInit(ami_mem, udata->nplist, udata->sensi_meth, fsxdot, sx, sdx);
    }

    int wrap_SetDenseJacFn() {
        return IDADlsSetDenseJacFn(ami_mem, fJ);
    }

    int wrap_SetSparseJacFn() {
        return IDASlsSetSparseJacFn(ami_mem, fJSparse);
    }

    int wrap_SetBandJacFn() {
        return IDADlsSetBandJacFn(ami_mem, fJBand);
    }

    int wrap_SetJacTimesVecFn() {
        return IDASpilsSetJacTimesVecFn(ami_mem, fJv);
    }

    int wrap_SetDenseJacFnB(int which) {
        return IDADlsSetDenseJacFnB(ami_mem, which, fJB);
    }

    int wrap_SetSparseJacFnB(int which) {
        return IDASlsSetSparseJacFnB(ami_mem, which, fJSparseB);
    }

    int wrap_SetBandJacFnB(int which) {
        return IDADlsSetBandJacFnB(ami_mem, which, fJBandB);
    }

    int wrap_SetJacTimesVecFnB(int which) {
        return IDASpilsSetJacTimesVecFnB(ami_mem, which, fJvB);
    }

    void *AMICreate(int lmm, int iter) {
        return IDACreate();
    }

    int AMISStolerances(double rtol,double atol) {
        return IDASStolerances(ami_mem,rtol,atol);
    }

    int AMISensEEtolerances() {
        return IDASensEEtolerances(mem);
    }

    int AMISetSensErrCon(bool error_corr) {
        return IDASetSensErrCon(ami_mem,error_corr);
    }

    int AMISetQuadErrConB(int which, bool flag) {
        return IDASetQuadErrConB(ami_mem,which,flag);
    }

    int AMIGetRootInfo(int *rootsfound) {
        return IDAGetRootInfo(ami_mem,rootsfound);
    }

    int AMISetErrHandlerFn() {
        return IDASetErrHandlerFn(ami_mem,wrap_ErrHandlerFn,NULL);
    }

    int AMISetUserData(void *user_data) {
        return IDASetUserData(ami_mem,user_data);
    }

    int AMISetUserDataB(int which, void *user_data) {
        return IDASetUserDataB(ami_mem, which, user_data);
    }

    int AMISetMaxNumSteps(long int mxsteps) {
        return IDASetMaxNumSteps(ami_mem,mxsteps);
    }

    int AMISetStabLimDet(int stldet) {
        return(0);
    }

    int AMISetStabLimDetB(int which, int stldet) {
        return(0);
    }

    int AMISetId(N_Vector id) {
        return IDASetId(ami_mem, id);
    }

    int AMISetSuppressAlg(bool flag) {
        return IDASetSuppressAlg(ami_mem, flag);
    }

    int AMIReInit(realtype t0, N_Vector yy0, N_Vector yp0) {
        return IDAReInit( ami_mem, t0, yy0, yp0);
    }

    int AMISensReInit(int ism, N_Vector *yS0, N_Vector *ypS0) {
        return IDASensReInit( ami_mem, ism, yS0, ypS0);
    }

    int AMISetSensParams(realtype *p, realtype *pbar, int *plist) {
        return IDASetSensParams(ami_mem, p, pbar, plist);
    }

    int AMIGetDky(realtype t, int k, N_Vector dky) {
        return IDAGetDky(ami_mem, t, k, dky);
    }

    int AMIGetSens(realtype *tret, N_Vector *yySout) {
        return IDAGetSens( ami_mem, tret, yySout);
    }

    int AMIRootInit(int nrtfn, IDARootFn ptr) {
        return IDARootInit( ami_mem, nrtfn, ptr);
    }

    void AMIFree() {
        IDAFree(ami_mem);
        ami_mem = NULL;
    }

    int AMIAdjInit(long int steps, int interp) {
        return IDAAdjInit(ami_mem, steps, interp);
    }

    int AMICreateB(int lmm, int iter, int *which) {
        return IDACreateB(ami_mem, which);
    }

    int AMIReInitB(int which, realtype tB0, N_Vector yyB0, N_Vector ypB0) {
        return IDAReInitB(ami_mem, which, tB0, yyB0, ypB0);
    }

    int AMISStolerancesB(int which, realtype relTolB, realtype absTolB) {
        return IDASStolerancesB(ami_mem, which, relTolB, absTolB);
    }

    int AMIQuadReInitB(int which, N_Vector yQB0) {
        return IDAQuadReInitB(ami_mem, which, yQB0);
    }

    int AMIQuadSStolerancesB(int which, realtype reltolQB, realtype abstolQB) {
        return IDAQuadSStolerancesB(ami_mem, which, reltolQB, abstolQB);
    }

    int AMISolve(realtype tout, N_Vector yret, N_Vector ypret, realtype *tret, int itask) {
        return IDASolve(ami_mem, tout, tret, yret, ypret, itask);
    }

    int AMISolveF(realtype tout, N_Vector yret, N_Vector ypret, realtype *tret, int itask, int *ncheckPtr) {
        return IDASolveF(ami_mem, tout, tret, yret, ypret, itask, ncheckPtr);
    }

    int AMISolveB(realtype tBout, int itaskB) {
        return IDASolveB(ami_mem, tBout, itaskB);
    }

    int AMISetMaxNumStepsB(int which, long int mxstepsB) {
        return IDASetMaxNumStepsB(ami_mem, which, mxstepsB);
    }

    int AMIGetB(int which, realtype *tret, N_Vector yy, N_Vector yp) {
        return IDAGetB(ami_mem, which, tret, yy, yp);
    }

    int AMIGetQuadB(int which, realtype *tret, N_Vector qB) {
        return IDAGetQuadB(ami_mem, which, tret, qB);
    }

    int AMIDense(int nx) {
        return IDADense(ami_mem, nx);
    }

    int AMIDenseB(int which, int nx) {
        return IDADenseB(ami_mem, which, nx);
    }

    int AMIBand(int nx, int ubw, int lbw) {
        return IDABand(ami_mem, nx, ubw, lbw);
    }

    int AMIBandB(int which, int nx, int ubw, int lbw) {
        return IDABandB(ami_mem, which, nx, ubw, lbw);
    }

    int AMIDiag() {
        return(-99);
    }

    int AMIDiagB(int which) {
        return(-99);
    }

    int AMISpgmr(int prectype, int maxl) {
        return IDASpgmr(ami_mem, maxl);
    }

    int AMISpgmrB(int which, int prectype, int maxl) {
        return IDASpgmrB(ami_mem, which, maxl);
    }

    int AMISpbcg(int prectype, int maxl) {
        return IDASpbcg(ami_mem, maxl);
    }

    int AMISpbcgB(int which, int prectype, int maxl) {
        return IDASpbcgB(ami_mem, which, maxl);
    }

    int AMISptfqmr(int prectype, int maxl) {
        return IDASptfqmr(ami_mem, maxl);
    }

    int AMISptfqmrB(int which, int prectype, int maxl) {
        return IDASptfqmrB(ami_mem, which, maxl);
    }

    int AMIKLU(int nx, int nnz, int sparsetype) {
        return IDAKLU(ami_mem, nx, nnz, sparsetype);
    }

    int AMIKLUSetOrdering(int ordering) {
        return IDAKLUSetOrdering(ami_mem, ordering);
    }

    int AMIKLUSetOrderingB(int which, int ordering) {
        return IDAKLUSetOrderingB(ami_mem, which, ordering);
    }

    int AMIKLUB(int which, int nx, int nnz, int sparsetype) {
        return IDAKLUB(ami_mem, which, nx, nnz, sparsetype);
    }

    int AMIGetNumSteps(void *ami_mem, long int *numsteps) {
        return IDAGetNumSteps(ami_mem,numsteps);
    }

    int AMIGetNumRhsEvals(void *ami_mem, long int *numrhsevals) {
        return IDAGetNumResEvals(ami_mem,numrhsevals);
    }

    int AMIGetNumErrTestFails(void *ami_mem, long int *numerrtestfails) {
        return IDAGetNumErrTestFails(ami_mem,numerrtestfails);
    }

    int AMIGetNumNonlinSolvConvFails(void *ami_mem, long int *numnonlinsolvconvfails) {
        return IDAGetNumNonlinSolvConvFails(ami_mem,numnonlinsolvconvfails);
    }

    int AMIGetLastOrder(int *order) {
        return IDAGetLastOrder(ami_mem,order);
    }

    void *AMIGetAdjBmem(void *ami_mem, int which) {
        return IDAGetAdjIDABmem(ami_mem,which);
    }

    int AMICalcIC(realtype tout1){
        return IDACalcIC(ami_mem,IDA_YA_YDP_INIT , tout1);
    }

    int AMICalcICB(int which, realtype tout1, N_Vector xB, N_Vector dxB) {
        return IDACalcICB(ami_mem, which, tout1, xB, dxB);
    }

    int AMISetStopTime(realtype tstop) {
        return IDASetStopTime(ami_mem, tstop);
    }

    static int resultFunction(realtype tt, N_Vector yy, N_Vector yp,
                          N_Vector rr, void *user_data) {
        return fxdot(tt, yy, yp, rr, user_data);
    }

    static int resultFunctionB(realtype tt,
                               N_Vector yy, N_Vector yp,
                               N_Vector yyB, N_Vector ypB,
                                           N_Vector rrB, void *user_dataB) {
        return fxBdot(tt, yy, yp, yyB, ypB, rrB, user_dataB);
    }


    ~IDASolver() {}
}

#endif /* idawrap_h */
