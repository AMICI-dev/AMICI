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
#include <include/amici_exception.h>
#include <include/tdata.h>
#include <include/udata.h>

namespace amici {

IDASolver::IDASolver() : Solver() {}

void IDASolver::init(N_Vector x, N_Vector dx, realtype t) {
    int status = IDAInit(ami_mem, residualFunction, RCONST(t), x, dx);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAInit");
    return;
}
void IDASolver::binit(int which, N_Vector xB, N_Vector dxB, realtype t) {
    int status = IDAInitB(ami_mem, which, residualFunctionB, RCONST(t), xB, dxB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAInitB");
    return;
}
void IDASolver::qbinit(int which, N_Vector qBdot) {
    int status = IDAQuadInitB(ami_mem, which, fqBdot, qBdot);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAQuadInitB");
    return;
}
void IDASolver::rootInit(int ne) {
    int status = IDARootInit(ami_mem, ne, rootFunction);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDARootInit");
    return;
}
void IDASolver::sensInit1(N_Vector *sx, N_Vector *sdx, const UserData *udata) {
    int status = IDASensInit(ami_mem, udata->nplist, udata->sensi_meth, fsxdot, sx,
                       sdx);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASensInit");
    return;
}
void IDASolver::setDenseJacFn() {  int status = IDADlsSetDenseJacFn(ami_mem, J); }

void IDASolver::setSparseJacFn() {
    int status = IDASlsSetSparseJacFn(ami_mem, fJSparse);
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDASlsSetSparseJacFn");
    return;
}
void IDASolver::setBandJacFn() {  int status = IDADlsSetBandJacFn(ami_mem, fJBand); }

void IDASolver::setJacTimesVecFn() {
    int status = IDASpilsSetJacTimesVecFn(ami_mem, fJv);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASpilsSetJacTimesVecFn");
    return;
}
void IDASolver::setDenseJacFnB(int which) {
    int status = IDADlsSetDenseJacFnB(ami_mem, which, fJB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDADlsSetDenseJacFnB");
    return;
}
void IDASolver::setSparseJacFnB(int which) {
    int status = IDASlsSetSparseJacFnB(ami_mem, which, fJSparseB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASlsSetSparseJacFnB");
    return;
}
void IDASolver::setBandJacFnB(int which) {
    int status = IDADlsSetBandJacFnB(ami_mem, which, fJBandB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDADlsSetBandJacFnB");
    return;
}
void IDASolver::setJacTimesVecFnB(int which) {
    int status = IDASpilsSetJacTimesVecFnB(ami_mem, which, fJvB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASpilsSetJacTimesVecFnB");
    return;
}
void *IDASolver::AMICreate(int lmm, int iter) {
    return IDACreate();
}

void IDASolver::AMISStolerances(double rtol, double atol) {
    int status = IDASStolerances(ami_mem, rtol, atol);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASStolerances");
    return;
}
void IDASolver::AMISensEEtolerances() {
    int status = IDASensEEtolerances(ami_mem);
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDASensEEtolerances");
    return;
}

void IDASolver::AMISetSensErrCon(bool error_corr) {
    int status = IDASetSensErrCon(ami_mem, error_corr);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetSensErrCon");
    return;
}
void IDASolver::AMISetQuadErrConB(int which, bool flag) {
    int status = IDASetQuadErrConB(ami_mem, which, flag);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetQuadErrConB");
    return;
}
void IDASolver::AMIGetRootInfo(int *rootsfound) {
    int status = IDAGetRootInfo(ami_mem, rootsfound);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetRootInfo");
    return;
}
void IDASolver::AMISetErrHandlerFn() {
    int status = IDASetErrHandlerFn(ami_mem, wrapErrHandlerFn, NULL);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetErrHandlerFn");
    return;
}
void IDASolver::AMISetUserData(void *user_data) {
    int status = IDASetUserData(ami_mem, user_data);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetUserData");
    return;
}
void IDASolver::AMISetUserDataB(int which, void *user_data) {
    int status = IDASetUserDataB(ami_mem, which, user_data);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetUserDataB");
    return;
}
void IDASolver::AMISetMaxNumSteps(long mxsteps) {
    int status = IDASetMaxNumSteps(ami_mem, mxsteps);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetMaxNumSteps");
    return;
}
void IDASolver::AMISetStabLimDet(int stldet) {
    return;
}

void IDASolver::AMISetStabLimDetB(int which, int stldet) {
    return;
}

void IDASolver::AMISetId(Model *model) {
    if (!model->idlist)
        throw AmiException("Model was not properly set up, missing definition of idlist");

    N_Vector id = N_VNew_Serial(model->nx);
    memcpy(NV_CONTENT_S(id)->data, model->idlist, model->nx * sizeof(realtype));

    int status = IDASetId(ami_mem, id);
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDASetMaxNumSteps");

    N_VDestroy_Serial(id);

    return;
}

void IDASolver::AMISetSuppressAlg(bool flag) {
    int status = IDASetSuppressAlg(ami_mem, flag);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetSuppressAlg");
    return;
}
void IDASolver::AMIReInit(realtype t0, N_Vector yy0, N_Vector yp0) {
    int status = IDAReInit(ami_mem, t0, yy0, yp0);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAReInit");
    return;
}
void IDASolver::AMISensReInit(int ism, N_Vector *yS0, N_Vector *ypS0) {
    int status = IDASensReInit(ami_mem, ism, yS0, ypS0);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASensReInit");
    return;
}
void IDASolver::AMISetSensParams(realtype *p, realtype *pbar, int *plist) {
    int status = IDASetSensParams(ami_mem, p, pbar, plist);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetSensParams");
    return;
}
void IDASolver::AMIGetDky(realtype t, int k, N_Vector dky) {
    int status = IDAGetDky(ami_mem, t, k, dky);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetDky");
    return;
}
void IDASolver::AMIGetSens(realtype *tret, N_Vector *yySout) {
    int status = IDAGetSens(ami_mem, tret, yySout);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetSens");
    return;
}
void IDASolver::AMIFree() {
    IDAFree(&ami_mem);
    ami_mem = NULL;
}

void IDASolver::AMIAdjInit(long steps, int interp) {
    int status = IDAAdjInit(ami_mem, steps, interp);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAAdjInit");
    return;
}
void IDASolver::AMICreateB(int lmm, int iter, int *which) {
    int status = IDACreateB(ami_mem, which);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDACreateB");
    return;
}
void IDASolver::AMIReInitB(int which, realtype tB0, N_Vector yyB0,
                          N_Vector ypB0) {
    int status = IDAReInitB(ami_mem, which, tB0, yyB0, ypB0);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAReInitB");
    return;
}
void IDASolver::AMISStolerancesB(int which, realtype relTolB, realtype absTolB) {
    int status = IDASStolerancesB(ami_mem, which, relTolB, absTolB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASStolerancesB");
    return;
}
void IDASolver::AMIQuadReInitB(int which, N_Vector yQB0) {
    int status = IDAQuadReInitB(ami_mem, which, yQB0);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAQuadReInitB");
    return;
}
void IDASolver::AMIQuadSStolerancesB(int which, realtype reltolQB,
                                    realtype abstolQB) {
    int status = IDAQuadSStolerancesB(ami_mem, which, reltolQB, abstolQB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAQuadSStolerancesB");
    return;
}
int IDASolver::AMISolve(realtype tout, N_Vector yret, N_Vector ypret,
                        realtype *tret, int itask) {
    return IDASolve(ami_mem, tout, tret, yret, ypret, itask);
}
int IDASolver::AMISolveF(realtype tout, N_Vector yret, N_Vector ypret,
                         realtype *tret, int itask, int *ncheckPtr) {
    return IDASolveF(ami_mem, tout, tret, yret, ypret, itask, ncheckPtr);
}
void IDASolver::AMISolveB(realtype tBout, int itaskB) {
    int status = IDASolveB(ami_mem, tBout, itaskB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASolveB");
    return;
}
void IDASolver::AMISetMaxNumStepsB(int which, long mxstepsB) {
    int status = IDASetMaxNumStepsB(ami_mem, which, mxstepsB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetMaxNumStepsB");
    return;
}
void IDASolver::AMIGetB(int which, realtype *tret, N_Vector yy, N_Vector yp) {
    int status = IDAGetB(ami_mem, which, tret, yy, yp);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetB");
    return;
}
void IDASolver::AMIGetQuadB(int which, realtype *tret, N_Vector qB) {
    int status = IDAGetQuadB(ami_mem, which, tret, qB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetQuadB");
    return;
}
void IDASolver::AMIDense(int nx) {
    int status = IDADense(ami_mem, nx);
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDADense");
    return;
}

void IDASolver::AMIDenseB(int which, int nx) {
    int status = IDADenseB(ami_mem, which, nx);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDADenseB");
    return;
}
                                            
void IDASolver::AMIBand(int nx, int ubw, int lbw) {
    int status = IDABand(ami_mem, nx, ubw, lbw);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDABand");
    return;
}
void IDASolver::AMIBandB(int which, int nx, int ubw, int lbw) {
    int status = IDABandB(ami_mem, which, nx, ubw, lbw);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDABandB");
    return;
}
void IDASolver::AMIDiag() {
    throw AmiException("Diag Solver was not implemented for DAEs");
    return;
}

void IDASolver::AMIDiagB(int which) {
    throw AmiException("Diag Solver was not implemented for DAEs");
    return;
}

void IDASolver::AMISpgmr(int prectype, int maxl) {
    int status = IDASpgmr(ami_mem, maxl);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASpgmr");
    return;
}
void IDASolver::AMISpgmrB(int which, int prectype, int maxl) {
    int status = IDASpgmrB(ami_mem, which, maxl);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASpgmrB");
    return;
}
void IDASolver::AMISpbcg(int prectype, int maxl) {
    int status = IDASpbcg(ami_mem, maxl);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASpbcg");
    return;
}
void IDASolver::AMISpbcgB(int which, int prectype, int maxl) {
    int status = IDASpbcgB(ami_mem, which, maxl);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASpbcgB");
    return;
}
void IDASolver::AMISptfqmr(int prectype, int maxl) {
    int status = IDASptfqmr(ami_mem, maxl);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASptfqmr");
    return;
}
void IDASolver::AMISptfqmrB(int which, int prectype, int maxl) {
    int status = IDASptfqmrB(ami_mem, which, maxl);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASptfqmrB");
    return;
}
void IDASolver::AMIKLU(int nx, int nnz, int sparsetype) {
    int status = IDAKLU(ami_mem, nx, nnz, sparsetype);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAKLU");
    return;
}
void IDASolver::AMIKLUSetOrdering(int ordering) {
    int status = IDAKLUSetOrdering(ami_mem, ordering);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAKLUSetOrdering");
    return;
}
void IDASolver::AMIKLUSetOrderingB(int which, int ordering) {
    int status = IDAKLUSetOrderingB(ami_mem, which, ordering);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAKLUSetOrderingB");
    return;
}
void IDASolver::AMIKLUB(int which, int nx, int nnz, int sparsetype) {
    int status = IDAKLUB(ami_mem, which, nx, nnz, sparsetype);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAKLUB");
    return;
}
void IDASolver::AMIGetNumSteps(void *ami_mem, long *numsteps) {
    int status = IDAGetNumSteps(ami_mem, numsteps);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetNumSteps");
    return;
}
void IDASolver::AMIGetNumRhsEvals(void *ami_mem, long *numrhsevals) {
    int status = IDAGetNumResEvals(ami_mem, numrhsevals);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetNumResEvals");
    return;
}
void IDASolver::AMIGetNumErrTestFails(void *ami_mem, long *numerrtestfails) {
    int status = IDAGetNumErrTestFails(ami_mem, numerrtestfails);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetNumErrTestFails");
    return;
}
void IDASolver::AMIGetNumNonlinSolvConvFails(void *ami_mem,
                                            long *numnonlinsolvconvfails) {
    int status = IDAGetNumNonlinSolvConvFails(ami_mem, numnonlinsolvconvfails);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetNumNonlinSolvConvFails");
    return;
}
void IDASolver::AMIGetLastOrder(void *ami_mem, int *order) {
    int status = IDAGetLastOrder(ami_mem, order);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetLastOrder");
    return;
}
void *IDASolver::AMIGetAdjBmem(void *ami_mem, int which) {
    return IDAGetAdjIDABmem(ami_mem, which);
}
void IDASolver::AMICalcIC(realtype tout1) {
    int status = IDACalcIC(ami_mem, IDA_YA_YDP_INIT, tout1);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDACalcIC");
    return;
}
void IDASolver::AMICalcICB(int which, realtype tout1, N_Vector xB,
                          N_Vector dxB) {
    int status = IDACalcICB(ami_mem, which, tout1, xB, dxB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDACalcICB");
    return;
}
void IDASolver::AMISetStopTime(realtype tstop) {
    int status = IDASetStopTime(ami_mem, tstop);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetStopTime");
    return;
}
                                            
void IDASolver::turnOffRootFinding() {
    int status = IDARootInit(ami_mem, 0, NULL);
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDARootInit");
    return;
}

int IDASolver::residualFunction(realtype tt, N_Vector yy, N_Vector yp,
                              N_Vector rr, void *user_data) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fxdot(tt, yy, yp, rr, user_data);
}
int IDASolver::residualFunctionB(realtype tt, N_Vector yy, N_Vector yp,
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
    return tdata->model->fqBdot(t, x, dx, xB, dxB, qBdot, user_data);
}
                                            
int IDASolver::fsxdot(int Ns, realtype t, N_Vector x, N_Vector xdot,
                      N_Vector dx, N_Vector *sx, N_Vector *sdx, N_Vector *sxdot,
                      void *user_data, N_Vector tmp1, N_Vector tmp2,
                      N_Vector tmp3) {
    TempData *tdata = (TempData *)user_data;
    int status;
    for(int ip = 0; ip<tdata->udata->nplist; ip++) {
        status = tdata->model->fsxdot(Ns, t, x, dx, xdot, ip, sx[tdata->udata->plist[ip]],
                                      sdx[tdata->udata->plist[ip]], sxdot[tdata->udata->plist[ip]],
                                      user_data, tmp1, tmp2, tmp3);
        if(status != AMICI_SUCCESS)
            return status;
    }
    return status;
}

int IDASolver::fJSparse(realtype t, realtype cj, N_Vector x, N_Vector dx,
                        N_Vector xdot, SlsMat J, void *user_data, N_Vector tmp1,
                        N_Vector tmp2, N_Vector tmp3) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fJSparse(t, cj, x, dx, xdot, J, user_data, tmp1, tmp2, tmp3);
}
                                            
int IDASolver::fJBand(long int N, long int mupper, long int mlower, realtype t,
                      realtype cj, N_Vector x, N_Vector dx, N_Vector xdot,
                      DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2,
                      N_Vector tmp3) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fJBand(N, mupper, mlower, t, cj, x, dx, xdot, J, user_data,
                                tmp1, tmp2, tmp3);
}
                                            
int IDASolver::fJv(realtype t, N_Vector x, N_Vector dx, N_Vector xdot,
                   N_Vector v, N_Vector Jv, realtype cj, void *user_data,
                   N_Vector tmp1, N_Vector tmp2) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fJv(t, x, dx, xdot, v, Jv, cj, user_data, tmp1, tmp2);
}
                                            
int IDASolver::fJB(long int NeqBdot, realtype t, realtype cj, N_Vector x,
                   N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot,
                   DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B,
                   N_Vector tmp3B) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fJB(NeqBdot, t, cj, x, dx, xB, dxB, xBdot, JB, user_data, tmp1B,
                             tmp2B, tmp3B);
}
                                            
int IDASolver::fJSparseB(realtype t, realtype cj, N_Vector x, N_Vector dx,
                         N_Vector xB, N_Vector dxB, N_Vector xBdot, SlsMat JB,
                         void *user_data, N_Vector tmp1B, N_Vector tmp2B,
                         N_Vector tmp3B) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fJSparseB(t, cj, x, dx, xB, dxB, xBdot, JB, user_data, tmp1B, tmp2B,
                                   tmp3B);
}
                                            
int IDASolver::fJBandB(long int NeqBdot, long int mupper, long int mlower,
                       realtype t, realtype cj, N_Vector x, N_Vector dx,
                       N_Vector xB, N_Vector dxB, N_Vector xBdot, DlsMat JB,
                       void *user_data, N_Vector tmp1B, N_Vector tmp2B,
                       N_Vector tmp3B) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fJBandB(NeqBdot, mupper, mlower, t, cj, x, dx, xB, dxB, xBdot, JB,
                                 user_data, tmp1B, tmp2B, tmp3B);
}
                                            
int IDASolver::fJvB(realtype t, N_Vector x, N_Vector dx, N_Vector xB,
                    N_Vector dxB, N_Vector xBdot, N_Vector vB, N_Vector JvB,
                    realtype cj, void *user_data, N_Vector tmpB1,
                    N_Vector tmpB2) {
    TempData *tdata = (TempData *)user_data;
    return tdata->model->fJvB(t, x, dx, xB, dxB, xBdot, vB, JvB, cj, user_data, tmpB1, tmpB2);
}
                                            
IDASolver::~IDASolver() { AMIFree(); }

} // namespace amici
