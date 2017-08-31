#include <include/amici_solver_idas.h>

#include <idas/idas.h>
/*#include <idas/idas_lapack.h>*/
#include <idas/idas_band.h>
#include <idas/idas_bbdpre.h>
#include <idas/idas_dense.h>
/*#include <idas/idas_diag.h>*/
#include <idas/idas_klu.h>
#include <idas/idas_spbcgs.h>
#include <idas/idas_spgmr.h>
#include <idas/idas_sptfqmr.h>

#include <amd.h>
#include <btf.h>
#include <colamd.h>
#include <klu.h>

#include <cstring>
#include <include/amici.h>
#include <include/amici_model.h>
#include <include/tdata.h>
#include <include/udata.h>

IDASolver::IDASolver() : Solver() {}

int IDASolver::init(N_Vector x, N_Vector dx, realtype t) {
    return IDAInit(ami_mem, resultFunction, RCONST(t), x, dx);
}

int IDASolver::binit(int which, N_Vector xB, N_Vector dxB, realtype t) {
    return IDAInitB(ami_mem, which, resultFunctionB, RCONST(t), xB, dxB);
}

int IDASolver::qbinit(int which, N_Vector qBdot) {
    return IDAQuadInitB(ami_mem, which, fqBdot, qBdot);
}

int IDASolver::rootInit(int ne) {
    return IDARootInit(ami_mem, ne, rootFunction);
}

int IDASolver::sensInit1(N_Vector *sx, N_Vector *sdx, const UserData *udata) {
    return IDASensInit(ami_mem, udata->nplist, udata->sensi_meth, fsxdot, sx,
                       sdx);
}

int IDASolver::setDenseJacFn() { return IDADlsSetDenseJacFn(ami_mem, J); }

int IDASolver::setSparseJacFn() {
    return IDASlsSetSparseJacFn(ami_mem, fJSparse);
}

int IDASolver::setBandJacFn() { return IDADlsSetBandJacFn(ami_mem, fJBand); }

int IDASolver::setJacTimesVecFn() {
    return IDASpilsSetJacTimesVecFn(ami_mem, fJv);
}

int IDASolver::setDenseJacFnB(int which) {
    return IDADlsSetDenseJacFnB(ami_mem, which, fJB);
}

int IDASolver::setSparseJacFnB(int which) {
    return IDASlsSetSparseJacFnB(ami_mem, which, fJSparseB);
}

int IDASolver::setBandJacFnB(int which) {
    return IDADlsSetBandJacFnB(ami_mem, which, fJBandB);
}

int IDASolver::setJacTimesVecFnB(int which) {
    return IDASpilsSetJacTimesVecFnB(ami_mem, which, fJvB);
}

void *IDASolver::AMICreate(int lmm, int iter) { return IDACreate(); }

int IDASolver::AMISStolerances(double rtol, double atol) {
    return IDASStolerances(ami_mem, rtol, atol);
}

int IDASolver::AMISensEEtolerances() { return IDASensEEtolerances(ami_mem); }

int IDASolver::AMISetSensErrCon(bool error_corr) {
    return IDASetSensErrCon(ami_mem, error_corr);
}

int IDASolver::AMISetQuadErrConB(int which, bool flag) {
    return IDASetQuadErrConB(ami_mem, which, flag);
}

int IDASolver::AMIGetRootInfo(int *rootsfound) {
    return IDAGetRootInfo(ami_mem, rootsfound);
}

int IDASolver::AMISetErrHandlerFn() {
    return IDASetErrHandlerFn(ami_mem, wrapErrHandlerFn, NULL);
}

int IDASolver::AMISetUserData(void *user_data) {
    return IDASetUserData(ami_mem, user_data);
}

int IDASolver::AMISetUserDataB(int which, void *user_data) {
    return IDASetUserDataB(ami_mem, which, user_data);
}

int IDASolver::AMISetMaxNumSteps(long mxsteps) {
    return IDASetMaxNumSteps(ami_mem, mxsteps);
}

int IDASolver::AMISetStabLimDet(int stldet) { return (0); }

int IDASolver::AMISetStabLimDetB(int which, int stldet) { return (0); }

int IDASolver::AMISetId(Model *model) {
    if (!model->idlist)
        return AMICI_ERROR_SETUP;

    N_Vector id = N_VNew_Serial(model->nx);
    if (!id)
        return AMICI_ERROR_SETUP;

    memcpy(NV_CONTENT_S(id)->data, model->idlist, model->nx * sizeof(realtype));

    int status = IDASetId(ami_mem, id);

    N_VDestroy_Serial(id);

    return status;
}

int IDASolver::AMISetSuppressAlg(bool flag) {
    return IDASetSuppressAlg(ami_mem, flag);
}

int IDASolver::AMIReInit(realtype t0, N_Vector yy0, N_Vector yp0) {
    return IDAReInit(ami_mem, t0, yy0, yp0);
}

int IDASolver::AMISensReInit(int ism, N_Vector *yS0, N_Vector *ypS0) {
    return IDASensReInit(ami_mem, ism, yS0, ypS0);
}

int IDASolver::AMISetSensParams(realtype *p, realtype *pbar, int *plist) {
    return IDASetSensParams(ami_mem, p, pbar, plist);
}

int IDASolver::AMIGetDky(realtype t, int k, N_Vector dky) {
    return IDAGetDky(ami_mem, t, k, dky);
}

int IDASolver::AMIGetSens(realtype *tret, N_Vector *yySout) {
    return IDAGetSens(ami_mem, tret, yySout);
}

void IDASolver::AMIFree() {
    IDAFree(&ami_mem);
    ami_mem = NULL;
}

int IDASolver::AMIAdjInit(long steps, int interp) {
    return IDAAdjInit(ami_mem, steps, interp);
}

int IDASolver::AMICreateB(int lmm, int iter, int *which) {
    return IDACreateB(ami_mem, which);
}

int IDASolver::AMIReInitB(int which, realtype tB0, N_Vector yyB0,
                          N_Vector ypB0) {
    return IDAReInitB(ami_mem, which, tB0, yyB0, ypB0);
}

int IDASolver::AMISStolerancesB(int which, realtype relTolB, realtype absTolB) {
    return IDASStolerancesB(ami_mem, which, relTolB, absTolB);
}

int IDASolver::AMIQuadReInitB(int which, N_Vector yQB0) {
    return IDAQuadReInitB(ami_mem, which, yQB0);
}

int IDASolver::AMIQuadSStolerancesB(int which, realtype reltolQB,
                                    realtype abstolQB) {
    return IDAQuadSStolerancesB(ami_mem, which, reltolQB, abstolQB);
}

int IDASolver::AMISolve(realtype tout, N_Vector yret, N_Vector ypret,
                        realtype *tret, int itask) {
    return IDASolve(ami_mem, tout, tret, yret, ypret, itask);
}

int IDASolver::AMISolveF(realtype tout, N_Vector yret, N_Vector ypret,
                         realtype *tret, int itask, int *ncheckPtr) {
    return IDASolveF(ami_mem, tout, tret, yret, ypret, itask, ncheckPtr);
}

int IDASolver::AMISolveB(realtype tBout, int itaskB) {
    return IDASolveB(ami_mem, tBout, itaskB);
}

int IDASolver::AMISetMaxNumStepsB(int which, long mxstepsB) {
    return IDASetMaxNumStepsB(ami_mem, which, mxstepsB);
}

int IDASolver::AMIGetB(int which, realtype *tret, N_Vector yy, N_Vector yp) {
    return IDAGetB(ami_mem, which, tret, yy, yp);
}

int IDASolver::AMIGetQuadB(int which, realtype *tret, N_Vector qB) {
    return IDAGetQuadB(ami_mem, which, tret, qB);
}

int IDASolver::AMIDense(int nx) { return IDADense(ami_mem, nx); }

int IDASolver::AMIDenseB(int which, int nx) {
    return IDADenseB(ami_mem, which, nx);
}

int IDASolver::AMIBand(int nx, int ubw, int lbw) {
    return IDABand(ami_mem, nx, ubw, lbw);
}

int IDASolver::AMIBandB(int which, int nx, int ubw, int lbw) {
    return IDABandB(ami_mem, which, nx, ubw, lbw);
}

int IDASolver::AMIDiag() { return (-99); }

int IDASolver::AMIDiagB(int which) { return (-99); }

int IDASolver::AMISpgmr(int prectype, int maxl) {
    return IDASpgmr(ami_mem, maxl);
}

int IDASolver::AMISpgmrB(int which, int prectype, int maxl) {
    return IDASpgmrB(ami_mem, which, maxl);
}

int IDASolver::AMISpbcg(int prectype, int maxl) {
    return IDASpbcg(ami_mem, maxl);
}

int IDASolver::AMISpbcgB(int which, int prectype, int maxl) {
    return IDASpbcgB(ami_mem, which, maxl);
}

int IDASolver::AMISptfqmr(int prectype, int maxl) {
    return IDASptfqmr(ami_mem, maxl);
}

int IDASolver::AMISptfqmrB(int which, int prectype, int maxl) {
    return IDASptfqmrB(ami_mem, which, maxl);
}

int IDASolver::AMIKLU(int nx, int nnz, int sparsetype) {
    return IDAKLU(ami_mem, nx, nnz, sparsetype);
}

int IDASolver::AMIKLUSetOrdering(int ordering) {
    return IDAKLUSetOrdering(ami_mem, ordering);
}

int IDASolver::AMIKLUSetOrderingB(int which, int ordering) {
    return IDAKLUSetOrderingB(ami_mem, which, ordering);
}

int IDASolver::AMIKLUB(int which, int nx, int nnz, int sparsetype) {
    return IDAKLUB(ami_mem, which, nx, nnz, sparsetype);
}

int IDASolver::AMIGetNumSteps(void *ami_mem, long *numsteps) {
    return IDAGetNumSteps(ami_mem, numsteps);
}

int IDASolver::AMIGetNumRhsEvals(void *ami_mem, long *numrhsevals) {
    return IDAGetNumResEvals(ami_mem, numrhsevals);
}

int IDASolver::AMIGetNumErrTestFails(void *ami_mem, long *numerrtestfails) {
    return IDAGetNumErrTestFails(ami_mem, numerrtestfails);
}

int IDASolver::AMIGetNumNonlinSolvConvFails(void *ami_mem,
                                            long *numnonlinsolvconvfails) {
    return IDAGetNumNonlinSolvConvFails(ami_mem, numnonlinsolvconvfails);
}

int IDASolver::AMIGetLastOrder(void *ami_mem, int *order) {
    return IDAGetLastOrder(ami_mem, order);
}

void *IDASolver::AMIGetAdjBmem(void *ami_mem, int which) {
    return IDAGetAdjIDABmem(ami_mem, which);
}

int IDASolver::AMICalcIC(realtype tout1) {
    return IDACalcIC(ami_mem, IDA_YA_YDP_INIT, tout1);
}

int IDASolver::AMICalcICB(int which, realtype tout1, N_Vector xB,
                          N_Vector dxB) {
    return IDACalcICB(ami_mem, which, tout1, xB, dxB);
}

int IDASolver::AMISetStopTime(realtype tstop) {
    return IDASetStopTime(ami_mem, tstop);
}

int IDASolver::turnOffRootFinding() { return IDARootInit(ami_mem, 0, NULL); }

int IDASolver::resultFunction(realtype tt, N_Vector yy, N_Vector yp,
                              N_Vector rr, void *user_data) {
    TempData *tdata = (TempData *)user_data;

    return tdata->model->fxdot(tt, yy, yp, rr, user_data);
}

int IDASolver::resultFunctionB(realtype tt, N_Vector yy, N_Vector yp,
                               N_Vector yyB, N_Vector ypB, N_Vector rrB,
                               void *user_dataB) {
    TempData *tdata = (TempData *)user_dataB;

    return tdata->model->fxBdot(tt, yy, yp, yyB, ypB, rrB, user_dataB);
}

int IDASolver::rootFunction(realtype t, N_Vector y, N_Vector yp, realtype *gout,
                            void *user_data) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->froot(t, y, yp, gout, user_data);
}

int IDASolver::J(long N, realtype t, realtype c_j, N_Vector y, N_Vector yp,
                 N_Vector r, DlsMat Jac, void *user_data, N_Vector tmp1,
                 N_Vector tmp2, N_Vector tmp3) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fJ(N, t, c_j, y, r, yp, Jac, user_data, tmp1, tmp2,
                            tmp3);
}

int IDASolver::fqBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB,
                      N_Vector dxB, N_Vector qBdot, void *user_data) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fqBdot(t, x, xB, qBdot, user_data);
}

int IDASolver::fsxdot(int Ns, realtype t, N_Vector x, N_Vector xdot,
                      N_Vector dx, N_Vector *sx, N_Vector *sxdot, N_Vector *sdx,
                      void *user_data, N_Vector tmp1, N_Vector tmp2,
                      N_Vector tmp3) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fsxdot(Ns, t, x, xdot, 0, *sx, *sxdot, user_data, tmp1,
                                tmp2);
}

int IDASolver::fJSparse(realtype t, realtype cj, N_Vector x, N_Vector dx,
                        N_Vector xdot, SlsMat J, void *user_data, N_Vector tmp1,
                        N_Vector tmp2, N_Vector tmp3) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fJSparse(t, x, xdot, J, user_data, tmp1, tmp2, tmp3);
}

int IDASolver::fJBand(long int N, long int mupper, long int mlower, realtype t,
                      realtype cj, N_Vector x, N_Vector dx, N_Vector xdot,
                      DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2,
                      N_Vector tmp3) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fJBand(N, mupper, mlower, t, x, xdot, J, user_data,
                                tmp1, tmp2, tmp3);
}

int IDASolver::fJv(realtype t, N_Vector x, N_Vector dx, N_Vector xdot,
                   N_Vector v, N_Vector Jv, realtype cj, void *user_data,
                   N_Vector tmp1, N_Vector tmp2) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fJv(v, Jv, t, x, xdot, user_data, tmp1);
}

int IDASolver::fJB(long int NeqBdot, realtype t, realtype cj, N_Vector x,
                   N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot,
                   DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B,
                   N_Vector tmp3B) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fJB(NeqBdot, t, x, xB, xBdot, JB, user_data, tmp1B,
                             tmp2B, tmp3B);
}

int IDASolver::fJSparseB(realtype t, realtype cj, N_Vector x, N_Vector dx,
                         N_Vector xB, N_Vector dxB, N_Vector xBdot, SlsMat JB,
                         void *user_data, N_Vector tmp1B, N_Vector tmp2B,
                         N_Vector tmp3B) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fJSparseB(t, x, xB, xBdot, JB, user_data, tmp1B, tmp2B,
                                   tmp3B);
}

int IDASolver::fJBandB(long int NeqBdot, long int mupper, long int mlower,
                       realtype t, realtype cj, N_Vector x, N_Vector dx,
                       N_Vector xB, N_Vector dxB, N_Vector xBdot, DlsMat JB,
                       void *user_data, N_Vector tmp1B, N_Vector tmp2B,
                       N_Vector tmp3B) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fJBandB(NeqBdot, mupper, mlower, t, x, xB, xBdot, JB,
                                 user_data, tmp1B, tmp2B, tmp3B);
}

int IDASolver::fJvB(realtype t, N_Vector x, N_Vector dx, N_Vector xB,
                    N_Vector dxB, N_Vector xBdot, N_Vector vB, N_Vector JvB,
                    realtype cj, void *user_data, N_Vector tmpB1,
                    N_Vector tmpB2) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fJvB(vB, JvB, t, x, xB, xBdot, user_data, tmpB1);
}

IDASolver::~IDASolver() { AMIFree(); }
