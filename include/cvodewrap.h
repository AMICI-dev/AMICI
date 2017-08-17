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
#include <include/amici_solver.h>
#include <include/amici_model_functions.h>
#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

class CVodeSolver : public Solver
{
public:
    CVodeSolver() {

    }

    int wrap_init(N_Vector x, N_Vector dx, realtype t) {
       return CVodeInit(ami_mem, resultFunction, RCONST(t), x);
    }

    int wrap_binit(int which, N_Vector xB, N_Vector dxB, realtype t) {
        return CVodeInitB(ami_mem, which, resultFunctionB, RCONST(t), xB);
    }

    int wrap_qbinit(int which, N_Vector qBdot) {
        return CVodeQuadInitB(ami_mem, which, fqBdot, qBdot);
    }

    int wrap_RootInit(UserData *udata) {
        return CVodeRootInit(ami_mem, udata->ne, rootFunction);
    }

    int wrap_SensInit1(N_Vector *sx, N_Vector *sdx, UserData *udata) {
        return CVodeSensInit1(ami_mem, udata->nplist, udata->sensi_meth, fsxdot, sx);
    }

    int wrap_SetDenseJacFn() {
        return CVDlsSetDenseJacFn(ami_mem, J);
    }

    int wrap_SetSparseJacFn() {
        return CVSlsSetSparseJacFn(ami_mem, fJSparse);
    }

    int wrap_SetBandJacFn() {
        return CVDlsSetBandJacFn(ami_mem, fJBand);
    }

    int wrap_SetJacTimesVecFn() {
        return CVSpilsSetJacTimesVecFn(ami_mem, fJv);
    }

    int wrap_SetDenseJacFnB(int which) {
        return CVDlsSetDenseJacFnB(ami_mem, which, fJB);
    }

    int wrap_SetSparseJacFnB(int which) {
        return CVSlsSetSparseJacFnB(ami_mem, which, fJSparseB);
    }

    int wrap_SetBandJacFnB(int which) {
        return CVDlsSetBandJacFnB(ami_mem, which, fJBandB);
    }

    int wrap_SetJacTimesVecFnB(int which) {
        return CVSpilsSetJacTimesVecFnB(ami_mem, which, fJvB);
    }

    void *AMICreate(int lmm, int iter) {
        return CVodeCreate(lmm,iter);
    }

    int AMISStolerances(double rtol,double atol) {
        return CVodeSStolerances(ami_mem,rtol,atol);
    }

    int AMISensEEtolerances() {
        return CVodeSensEEtolerances(ami_mem);
    }

    int AMISetSensErrCon(bool error_corr) {
        return CVodeSetSensErrCon(ami_mem,error_corr);
    }

    int AMISetQuadErrConB(int which, bool flag) {
        return CVodeSetQuadErrConB(ami_mem,which,flag);
    }

    int AMIGetRootInfo(int *rootsfound) {
        return CVodeGetRootInfo(ami_mem,rootsfound);
    }

    int AMISetErrHandlerFn() {
        return CVodeSetErrHandlerFn(ami_mem,wrap_ErrHandlerFn,NULL);
    }

    int AMISetUserData(void *user_data) {
        return CVodeSetUserData(ami_mem,user_data);
    }

    int AMISetUserDataB(int which, void *user_data) {
        return CVodeSetUserDataB(ami_mem, which, user_data);
    }

    int AMISetMaxNumSteps(long int mxsteps) {
        return CVodeSetMaxNumSteps(ami_mem,mxsteps);
    }

    int AMISetStabLimDet(int stldet) {
        return CVodeSetStabLimDet(ami_mem,stldet);
    }

    int AMISetStabLimDetB(int which, int stldet) {
        return CVodeSetStabLimDetB(ami_mem,which,stldet);
    }

    int AMISetId(N_Vector id) {
        return(0);
    }

    int AMISetSuppressAlg(bool flag) {
        return(0);
    }

    int AMIReInit(realtype t0, N_Vector yy0, N_Vector yp0) {
        return CVodeReInit( ami_mem, t0, yy0);
    }

    int AMISensReInit(int ism, N_Vector *yS0, N_Vector *ypS0) {
        return CVodeSensReInit( ami_mem, ism, yS0);
    }

    int AMISetSensParams(realtype *p, realtype *pbar, int *plist) {
        return CVodeSetSensParams(ami_mem, p, pbar, plist);
    }

    int AMIGetDky(realtype t, int k, N_Vector dky) {
        return CVodeGetDky(ami_mem, t, k, dky);
    }

    int AMIGetSens(realtype *tret, N_Vector *yySout) {
        return CVodeGetSens( ami_mem, tret, yySout);
    }

    int AMIRootInit(int nrtfn, CVRootFn ptr) {
        return CVodeRootInit( ami_mem, nrtfn, ptr);
    }

    void AMIFree(void **mem) {
        CVodeFree(mem);
        mem = NULL;
    }

    int AMIAdjInit(long int steps, int interp) {
        return CVodeAdjInit(ami_mem, steps, interp);
    }

    int AMICreateB(int lmm, int iter, int *which) {
        return CVodeCreateB(ami_mem, lmm, iter, which);
    }

    int AMIReInitB(int which, realtype tB0, N_Vector yyB0, N_Vector ypB0) {
        return CVodeReInitB(ami_mem, which, tB0, yyB0);
    }

    int AMISStolerancesB(int which, realtype relTolB, realtype absTolB) {
        return CVodeSStolerancesB(ami_mem, which, relTolB, absTolB);
    }

    int AMIQuadReInitB(int which, N_Vector yQB0) {
        return CVodeQuadReInitB(ami_mem, which, yQB0);
    }

    int AMIQuadSStolerancesB(int which, realtype reltolQB, realtype abstolQB) {
        return CVodeQuadSStolerancesB(ami_mem, which, reltolQB, abstolQB);
    }

    int AMISolve(realtype tout, N_Vector yret, N_Vector ypret, realtype *tret, int itask) {
        return CVode(ami_mem, tout, yret, tret, itask);
    }

    int AMISolveF(realtype tout, N_Vector yret, N_Vector ypret, realtype *tret, int itask, int *ncheckPtr) {
        return CVodeF(ami_mem, tout, yret, tret, itask, ncheckPtr);
    }

    int AMISolveB(realtype tBout, int itaskB) {
        return CVodeB(ami_mem, tBout, itaskB);
    }

    int AMISetMaxNumStepsB(int which, long int mxstepsB) {
        return CVodeSetMaxNumStepsB(ami_mem, which, mxstepsB);
    }

    int AMIGetB(int which, realtype *tret, N_Vector yy, N_Vector yp) {
        return CVodeGetB(ami_mem, which, tret, yy);
    }

    int AMIGetQuadB(int which, realtype *tret, N_Vector qB) {
        return CVodeGetQuadB(ami_mem, which, tret, qB);
    }

    int AMIDense(int nx) {
        return CVDense(ami_mem, nx);
    }

    int AMIDenseB(int which, int nx) {
        return CVDenseB(ami_mem, which, nx);
    }

    int AMIBand(int nx, int ubw, int lbw) {
        return CVBand(ami_mem, nx, ubw, lbw);
    }

    int AMIBandB(int which, int nx, int ubw, int lbw) {
        return CVBandB(ami_mem, which, nx, ubw, lbw);
    }

    int AMIDiag() {
        return CVDiag(ami_mem);
    }

    int AMIDiagB(int which) {
        return CVDiagB(ami_mem,which);
    }

    int AMISpgmr(int prectype, int maxl) {
        return CVSpgmr(ami_mem, prectype, maxl);
    }

    int AMISpgmrB(int which, int prectype, int maxl) {
        return CVSpgmrB(ami_mem, which, prectype, maxl);
    }

    int AMISpbcg(int prectype, int maxl) {
        return CVSpbcg(ami_mem, prectype, maxl);
    }

    int AMISpbcgB(int which, int prectype, int maxl) {
        return CVSpbcgB(ami_mem, which, prectype, maxl);
    }

    int AMISptfqmr(int prectype, int maxl) {
        return CVSptfqmr(ami_mem, prectype, maxl);
    }

    int AMISptfqmrB(int which, int prectype, int maxl) {
        return CVSptfqmrB(ami_mem, which, prectype, maxl);
    }

    int AMIKLU(int nx, int nnz, int sparsetype) {
        return CVKLU(ami_mem, nx, nnz, sparsetype);
    }

    int AMIKLUSetOrdering(int ordering) {
        return CVKLUSetOrdering(ami_mem, ordering);
    }

    int AMIKLUSetOrderingB(int which, int ordering) {
        return CVKLUSetOrderingB(ami_mem, which, ordering);
    }

    int AMIKLUB(int which, int nx, int nnz, int sparsetype) {
        return CVKLUB(ami_mem, which, nx, nnz, sparsetype);
    }

    int AMIGetNumSteps(void *ami_mem, long int *numsteps) {
        return CVodeGetNumSteps(ami_mem,numsteps);
    }

    int AMIGetNumRhsEvals(void *ami_mem, long int *numrhsevals) {
        return CVodeGetNumRhsEvals(ami_mem,numrhsevals);
    }

    int AMIGetNumErrTestFails(void *ami_mem, long int *numerrtestfails) {
        return CVodeGetNumErrTestFails(ami_mem,numerrtestfails);
    }

    int AMIGetNumNonlinSolvConvFails(void *ami_mem, long int *numnonlinsolvconvfails) {
        return CVodeGetNumNonlinSolvConvFails(ami_mem,numnonlinsolvconvfails);
    }

    int AMIGetLastOrder(void *ami_ami_mem, int *order) {
        return CVodeGetLastOrder(ami_mem,order);
    }

    void *AMIGetAdjBmem(void *ami_mem, int which) {
        return CVodeGetAdjCVodeBmem(ami_mem,which);
    }

    int AMICalcIC(realtype tout1) {
        return(0);
    }

    int AMICalcICB(int which, realtype tout1, N_Vector xB, N_Vector dxB) {
        return(0);
    }

    int AMISetStopTime(realtype tstop) {
        return CVodeSetStopTime(ami_mem, tstop);
    }

    static int resultFunction(realtype t, N_Vector y,
                   N_Vector ydot, void *user_data) {
        return fxdot(t, y, NULL, ydot, user_data);
    }

    static int resultFunctionB(realtype t, N_Vector y,
                                   N_Vector yB, N_Vector yBdot,
                                   void *user_dataB) {
        return fxBdot(t, y, yB, yBdot, user_dataB);
    }

    static int rootFunction(realtype t, N_Vector x, realtype *root, void *user_data) {
        return froot(t, x, NULL, root, user_data);
    }

    static int J(long int N, realtype t, N_Vector x, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
        return fJ(N, t, 0.0, x, NULL, xdot, J, user_data, tmp1, tmp2, tmp3);
    }

    ~CVodeSolver() {

    }
};

#endif /* CVodewrap_h */
