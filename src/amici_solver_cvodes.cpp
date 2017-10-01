#include "include/amici_solver_cvodes.h"

#include <cvodes/cvodes.h>
/*#include <cvodes/cvodes_lapack.h>*/
#include <cvodes/cvodes_band.h>
#include <cvodes/cvodes_bbdpre.h>
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_diag.h>
#include <cvodes/cvodes_klu.h>
#include <cvodes/cvodes_spbcgs.h>
#include <cvodes/cvodes_spgmr.h>
#include <cvodes/cvodes_sptfqmr.h>

#include <amd.h>
#include <btf.h>
#include <colamd.h>
#include <klu.h>

#include <include/amici.h>
#include <include/amici_model.h>
#include <include/amici_solver.h>
#include <include/tdata.h>
#include <include/udata.h>



int CVodeSolver::init(N_Vector x, N_Vector dx, realtype t) {
    return CVodeInit(ami_mem, resultFunction, RCONST(t), x);
}

int CVodeSolver::binit(int which, N_Vector xB, N_Vector dxB, realtype t) {
    return CVodeInitB(ami_mem, which, resultFunctionB, RCONST(t), xB);
}

int CVodeSolver::qbinit(int which, N_Vector qBdot) {
    return CVodeQuadInitB(ami_mem, which, fqBdot, qBdot);
}

int CVodeSolver::qbsinit(int which, N_Vector qBdot) {
    return CVodeQuadInitBS(ami_mem, which, fqBo2dot, qBdot);
}

int CVodeSolver::rootInit(int ne) {
    return CVodeRootInit(ami_mem, ne, rootFunction);
}

int CVodeSolver::sensInit1(N_Vector *sx, N_Vector *sdx, const UserData *udata) {
    return CVodeSensInit1(ami_mem, udata->nplist, udata->sensi_meth, fsxdot,
                          sx);
}

int CVodeSolver::setDenseJacFn() { return CVDlsSetDenseJacFn(ami_mem, J); }

int CVodeSolver::setSparseJacFn() {
    return CVSlsSetSparseJacFn(ami_mem, fJSparse);
}

int CVodeSolver::setBandJacFn() { return CVDlsSetBandJacFn(ami_mem, fJBand); }

int CVodeSolver::setJacTimesVecFn() {
    return CVSpilsSetJacTimesVecFn(ami_mem, fJv);
}

int CVodeSolver::setDenseJacFnB(int which) {
    return CVDlsSetDenseJacFnB(ami_mem, which, fJB);
}

int CVodeSolver::setSparseJacFnB(int which) {
    return CVSlsSetSparseJacFnB(ami_mem, which, fJSparseB);
}

int CVodeSolver::setBandJacFnB(int which) {
    return CVDlsSetBandJacFnB(ami_mem, which, fJBandB);
}

int CVodeSolver::setJacTimesVecFnB(int which) {
    return CVSpilsSetJacTimesVecFnB(ami_mem, which, fJvB);
}

void *CVodeSolver::AMICreate(int lmm, int iter) {
    return CVodeCreate(lmm, iter);
}

int CVodeSolver::AMISStolerances(double rtol, double atol) {
    return CVodeSStolerances(ami_mem, rtol, atol);
}

int CVodeSolver::AMISensEEtolerances() {
    return CVodeSensEEtolerances(ami_mem);
}

int CVodeSolver::AMISetSensErrCon(bool error_corr) {
    return CVodeSetSensErrCon(ami_mem, error_corr);
}

int CVodeSolver::AMISetQuadErrConB(int which, bool flag) {
    return CVodeSetQuadErrConB(ami_mem, which, flag);
}

int CVodeSolver::AMIGetRootInfo(int *rootsfound) {
    return CVodeGetRootInfo(ami_mem, rootsfound);
}

int CVodeSolver::AMISetErrHandlerFn() {
    return CVodeSetErrHandlerFn(ami_mem, wrapErrHandlerFn, NULL);
}

int CVodeSolver::AMISetUserData(void *user_data) {
    return CVodeSetUserData(ami_mem, user_data);
}

int CVodeSolver::AMISetUserDataB(int which, void *user_data) {
    return CVodeSetUserDataB(ami_mem, which, user_data);
}

int CVodeSolver::AMISetMaxNumSteps(long mxsteps) {
    return CVodeSetMaxNumSteps(ami_mem, mxsteps);
}

int CVodeSolver::AMISetStabLimDet(int stldet) {
    return CVodeSetStabLimDet(ami_mem, stldet);
}

int CVodeSolver::AMISetStabLimDetB(int which, int stldet) {
    return CVodeSetStabLimDetB(ami_mem, which, stldet);
}

int CVodeSolver::AMISetId(Model *model) { return (0); }

int CVodeSolver::AMISetSuppressAlg(bool flag) { return (0); }

int CVodeSolver::AMIReInit(realtype t0, N_Vector yy0, N_Vector yp0) {
    return CVodeReInit(ami_mem, t0, yy0);
}

int CVodeSolver::AMISensReInit(int ism, N_Vector *yS0, N_Vector *ypS0) {
    return CVodeSensReInit(ami_mem, ism, yS0);
}

int CVodeSolver::AMISetSensParams(realtype *p, realtype *pbar, int *plist) {
    return CVodeSetSensParams(ami_mem, p, pbar, plist);
}

int CVodeSolver::AMIGetDky(realtype t, int k, N_Vector dky) {
    return CVodeGetDky(ami_mem, t, k, dky);
}

int CVodeSolver::AMIGetSens(realtype *tret, N_Vector *yySout) {
    return CVodeGetSens(ami_mem, tret, yySout);
}

// int CVodeSolver::AMIRootInit(int nrtfn, RootFn ptr) {
//    return CVodeRootInit( ami_mem, nrtfn, ptr);
//}

void CVodeSolver::AMIFree() {
    CVodeFree(&ami_mem);
    ami_mem = NULL;
}

int CVodeSolver::AMIAdjInit(long steps, int interp) {
    return CVodeAdjInit(ami_mem, steps, interp);
}

int CVodeSolver::AMICreateB(int lmm, int iter, int *which) {
    return CVodeCreateB(ami_mem, lmm, iter, which);
}

int CVodeSolver::AMIReInitB(int which, realtype tB0, N_Vector yyB0,
                            N_Vector ypB0) {
    return CVodeReInitB(ami_mem, which, tB0, yyB0);
}

int CVodeSolver::AMISStolerancesB(int which, realtype relTolB,
                                  realtype absTolB) {
    return CVodeSStolerancesB(ami_mem, which, relTolB, absTolB);
}

int CVodeSolver::AMIQuadReInitB(int which, N_Vector yQB0) {
    return CVodeQuadReInitB(ami_mem, which, yQB0);
}

int CVodeSolver::AMIQuadSStolerancesB(int which, realtype reltolQB,
                                      realtype abstolQB) {
    return CVodeQuadSStolerancesB(ami_mem, which, reltolQB, abstolQB);
}

int CVodeSolver::AMISolve(realtype tout, N_Vector yret, N_Vector ypret,
                          realtype *tret, int itask) {
    return CVode(ami_mem, tout, yret, tret, itask);
}

int CVodeSolver::AMISolveF(realtype tout, N_Vector yret, N_Vector ypret,
                           realtype *tret, int itask, int *ncheckPtr) {
    return CVodeF(ami_mem, tout, yret, tret, itask, ncheckPtr);
}

int CVodeSolver::AMISolveB(realtype tBout, int itaskB) {
    return CVodeB(ami_mem, tBout, itaskB);
}

int CVodeSolver::AMISetMaxNumStepsB(int which, long mxstepsB) {
    return CVodeSetMaxNumStepsB(ami_mem, which, mxstepsB);
}

int CVodeSolver::AMIGetB(int which, realtype *tret, N_Vector yy, N_Vector yp) {
    return CVodeGetB(ami_mem, which, tret, yy);
}

int CVodeSolver::AMIGetQuadB(int which, realtype *tret, N_Vector qB) {
    return CVodeGetQuadB(ami_mem, which, tret, qB);
}

int CVodeSolver::AMIDense(int nx) { return CVDense(ami_mem, nx); }

int CVodeSolver::AMIDenseB(int which, int nx) {
    return CVDenseB(ami_mem, which, nx);
}

int CVodeSolver::AMIBand(int nx, int ubw, int lbw) {
    return CVBand(ami_mem, nx, ubw, lbw);
}

int CVodeSolver::AMIBandB(int which, int nx, int ubw, int lbw) {
    return CVBandB(ami_mem, which, nx, ubw, lbw);
}

int CVodeSolver::AMIDiag() { return CVDiag(ami_mem); }

int CVodeSolver::AMIDiagB(int which) { return CVDiagB(ami_mem, which); }

int CVodeSolver::AMISpgmr(int prectype, int maxl) {
    return CVSpgmr(ami_mem, prectype, maxl);
}

int CVodeSolver::AMISpgmrB(int which, int prectype, int maxl) {
    return CVSpgmrB(ami_mem, which, prectype, maxl);
}

int CVodeSolver::AMISpbcg(int prectype, int maxl) {
    return CVSpbcg(ami_mem, prectype, maxl);
}

int CVodeSolver::AMISpbcgB(int which, int prectype, int maxl) {
    return CVSpbcgB(ami_mem, which, prectype, maxl);
}

int CVodeSolver::AMISptfqmr(int prectype, int maxl) {
    return CVSptfqmr(ami_mem, prectype, maxl);
}

int CVodeSolver::AMISptfqmrB(int which, int prectype, int maxl) {
    return CVSptfqmrB(ami_mem, which, prectype, maxl);
}

int CVodeSolver::AMIKLU(int nx, int nnz, int sparsetype) {
    return CVKLU(ami_mem, nx, nnz, sparsetype);
}

int CVodeSolver::AMIKLUSetOrdering(int ordering) {
    return CVKLUSetOrdering(ami_mem, ordering);
}

int CVodeSolver::AMIKLUSetOrderingB(int which, int ordering) {
    return CVKLUSetOrderingB(ami_mem, which, ordering);
}

int CVodeSolver::AMIKLUB(int which, int nx, int nnz, int sparsetype) {
    return CVKLUB(ami_mem, which, nx, nnz, sparsetype);
}

int CVodeSolver::AMIGetNumSteps(void *ami_mem, long *numsteps) {
    return CVodeGetNumSteps(ami_mem, numsteps);
}

int CVodeSolver::AMIGetNumRhsEvals(void *ami_mem, long *numrhsevals) {
    return CVodeGetNumRhsEvals(ami_mem, numrhsevals);
}

int CVodeSolver::AMIGetNumErrTestFails(void *ami_mem, long *numerrtestfails) {
    return CVodeGetNumErrTestFails(ami_mem, numerrtestfails);
}

int CVodeSolver::AMIGetNumNonlinSolvConvFails(void *ami_mem,
                                              long *numnonlinsolvconvfails) {
    return CVodeGetNumNonlinSolvConvFails(ami_mem, numnonlinsolvconvfails);
}

int CVodeSolver::AMIGetLastOrder(void *ami_ami_mem, int *order) {
    return CVodeGetLastOrder(ami_mem, order);
}

void *CVodeSolver::AMIGetAdjBmem(void *ami_mem, int which) {
    return CVodeGetAdjCVodeBmem(ami_mem, which);
}

int CVodeSolver::AMICalcIC(realtype tout1) { return (0); }

int CVodeSolver::AMICalcICB(int which, realtype tout1, N_Vector xB,
                            N_Vector dxB) {
    return (0);
}

int CVodeSolver::AMISetStopTime(realtype tstop) {
    return CVodeSetStopTime(ami_mem, tstop);
}

int CVodeSolver::turnOffRootFinding() {
    return CVodeRootInit(ami_mem, 0, NULL);
}

int CVodeSolver::resultFunction(realtype t, N_Vector y, N_Vector ydot,
                                void *user_data) {
    TempData *tdata = (TempData *)user_data;

    return tdata->model->fxdot(t, y, NULL, ydot, user_data);
}

int CVodeSolver::resultFunctionB(realtype t, N_Vector y, N_Vector yB,
                                 N_Vector yBdot, void *user_dataB) {
    TempData *tdata = (TempData *)user_dataB;

    return tdata->model->fxBdot(t, y, NULL, yB, NULL, yBdot, user_dataB);
}

int CVodeSolver::rootFunction(realtype t, N_Vector x, realtype *root,
                              void *user_data) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->froot(t, x, NULL, root, user_data);
}

int CVodeSolver::J(long N, realtype t, N_Vector x, N_Vector xdot, DlsMat J,
                   void *user_data, N_Vector tmp1, N_Vector tmp2,
                   N_Vector tmp3) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fJ(N, t, 0.0, x, NULL, xdot, J, user_data, tmp1, tmp2,
                            tmp3);
}

int CVodeSolver::fqBdot(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot,
                        void *user_data) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fqBdot(t, x, xB, qBdot, user_data);
}

int CVodeSolver::fqBo2dot(realtype t, N_Vector x, N_Vector *sx, N_Vector xB, N_Vector qBdot,
                        void *user_data) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fqBo2dot(t, x, sx, xB, qBdot, user_data);
}

int CVodeSolver::fsxdot(int Ns, realtype t, N_Vector x, N_Vector xdot, int ip,
                        N_Vector sx, N_Vector sxdot, void *user_data,
                        N_Vector tmp1, N_Vector tmp2) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fsxdot(Ns, t, x, xdot, ip, sx, sxdot, user_data, tmp1,
                                tmp2);
}

int CVodeSolver::fJSparse(realtype t, N_Vector x, N_Vector xdot, SlsMat J,
                          void *user_data, N_Vector tmp1, N_Vector tmp2,
                          N_Vector tmp3) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fJSparse(t, x, xdot, J, user_data, tmp1, tmp2, tmp3);
}

int CVodeSolver::fJBand(long N, long mupper, long mlower, realtype t,
                        N_Vector x, N_Vector xdot, DlsMat J, void *user_data,
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fJBand(N, mupper, mlower, t, x, xdot, J, user_data,
                                tmp1, tmp2, tmp3);
}

int CVodeSolver::fJv(N_Vector v, N_Vector Jv, realtype t, N_Vector x,
                     N_Vector xdot, void *user_data, N_Vector tmp) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fJv(v, Jv, t, x, xdot, user_data, tmp);
}

int CVodeSolver::fJB(long int NeqBdot, realtype t, N_Vector x, N_Vector xB,
                     N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B,
                     N_Vector tmp2B, N_Vector tmp3B) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fJB(NeqBdot, t, x, xB, xBdot, JB, user_data, tmp1B,
                             tmp2B, tmp3B);
}

int CVodeSolver::fJSparseB(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot,
                           SlsMat JB, void *user_data, N_Vector tmp1B,
                           N_Vector tmp2B, N_Vector tmp3B) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fJSparseB(t, x, xB, xBdot, JB, user_data, tmp1B, tmp2B,
                                   tmp3B);
}

int CVodeSolver::fJBandB(long int NeqBdot, long int mupper, long int mlower,
                         realtype t, N_Vector x, N_Vector xB, N_Vector xBdot,
                         DlsMat JB, void *user_data, N_Vector tmp1B,
                         N_Vector tmp2B, N_Vector tmp3B) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fJBandB(NeqBdot, mupper, mlower, t, x, xB, xBdot, JB,
                                 user_data, tmp1B, tmp2B, tmp3B);
}

int CVodeSolver::fJvB(N_Vector vB, N_Vector JvB, realtype t, N_Vector x,
                      N_Vector xB, N_Vector xBdot, void *user_data,
                      N_Vector tmpB) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fJvB(vB, JvB, t, x, xB, xBdot, user_data, tmpB);
}

CVodeSolver::~CVodeSolver() { AMIFree(); }
