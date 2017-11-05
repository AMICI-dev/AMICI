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

#include <include/amici_model_ode.h>
#include <include/amici.h>
#include <include/amici_solver_cvodes.h>
#include <include/amici_exception.h>
#include <include/udata.h>

/**
 * @ brief extract information from a property of a matlab class (matrix)
 * @ param OPTION name of the property
 */

namespace amici {

void CVodeSolver::init(N_Vector x, N_Vector dx, realtype t) {
    int status = CVodeInit(ami_mem, model->fxdot, RCONST(t), x);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeInit");
}

void CVodeSolver::binit(int which, N_Vector xB, N_Vector dxB, realtype t) {
    int status = CVodeInitB(ami_mem, which, model->fxBdot, RCONST(t), xB);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeInitB");
}

void CVodeSolver::qbinit(int which, N_Vector qBdot) {
    int status = CVodeQuadInitB(ami_mem, which, model->fqBdot, qBdot);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeQuadInitB");
}

void CVodeSolver::rootInit(int ne) {
    int status = CVodeRootInit(ami_mem, ne, model->froot);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeRootInit");
}

void CVodeSolver::sensInit1(N_Vector *sx, N_Vector *sdx, const UserData *udata) {
    int status = CVodeSensInit1(ami_mem, udata->nplist(), udata->sensi_meth, model->fsxdot,
                          sx);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSensInit1");
}

void CVodeSolver::setDenseJacFn() {
    int status = CVDlsSetDenseJacFn(ami_mem, model->J);
    if(status != CV_SUCCESS)
        throw CvodeException(status,"CVDlsSetDenseJacFn");
}

void CVodeSolver::setSparseJacFn() {
    int status = CVSlsSetSparseJacFn(ami_mem, model->fJSparse);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVSlsSetSparseJacFn");
}

void CVodeSolver::setBandJacFn() {
    int status = CVDlsSetBandJacFn(ami_mem, model->fJBand);
    if(status != CV_SUCCESS)
        throw CvodeException(status,"CVDlsSetBandJacFn");
}

void CVodeSolver::setJacTimesVecFn() {
    int status = CVSpilsSetJacTimesVecFn(ami_mem, model->fJv);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVSpilsSetJacTimesVecFn");
}

void CVodeSolver::setDenseJacFnB(int which) {
    int status = CVDlsSetDenseJacFnB(ami_mem, which, model->fJB);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVDlsSetDenseJacFnB");
}

void CVodeSolver::setSparseJacFnB(int which) {
    int status = CVSlsSetSparseJacFnB(ami_mem, which, model->fJSparseB);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVSlsSetSparseJacFnB");
}

void CVodeSolver::setBandJacFnB(int which) {
    int status = CVDlsSetBandJacFnB(ami_mem, which, model->fJBandB);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVDlsSetBandJacFnB");
}

void CVodeSolver::setJacTimesVecFnB(int which) {
    int status = CVSpilsSetJacTimesVecFnB(ami_mem, which, model->fJvB);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVSpilsSetJacTimesVecFnB");
}

void *CVodeSolver::AMICreate(int lmm, int iter) {
    return CVodeCreate(lmm, iter);
}

void CVodeSolver::AMISStolerances(double rtol, double atol) {
    int status = CVodeSStolerances(ami_mem, rtol, atol);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSStolerances");
}

void CVodeSolver::AMISensEEtolerances() {
    int status = CVodeSensEEtolerances(ami_mem);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSensEEtolerances");
}

void CVodeSolver::AMISetSensErrCon(bool error_corr) {
    int status = CVodeSetSensErrCon(ami_mem, error_corr);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetSensErrCon");
}

void CVodeSolver::AMISetQuadErrConB(int which, bool flag) {
    int status = CVodeSetQuadErrConB(ami_mem, which, flag);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetQuadErrConB");
}

void CVodeSolver::AMIGetRootInfo(int *rootsfound) {
    int status = CVodeGetRootInfo(ami_mem, rootsfound);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeGetRootInfo");
}

void CVodeSolver::AMISetErrHandlerFn() {
    int status = CVodeSetErrHandlerFn(ami_mem, wrapErrHandlerFn, NULL);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetErrHandlerFn");
}

void CVodeSolver::AMISetUserData(void *user_data) {
    int status = CVodeSetUserData(ami_mem, user_data);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetUserData");
}

void CVodeSolver::AMISetUserDataB(int which, void *user_data) {
    int status = CVodeSetUserDataB(ami_mem, which, user_data);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetUserDataB");
}

void CVodeSolver::AMISetMaxNumSteps(long mxsteps) {
    int status = CVodeSetMaxNumSteps(ami_mem, mxsteps);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetMaxNumSteps");
}

void CVodeSolver::AMISetStabLimDet(int stldet) {
    int status = CVodeSetStabLimDet(ami_mem, stldet);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetStabLimDet");
}

void CVodeSolver::AMISetStabLimDetB(int which, int stldet) {
    int status = CVodeSetStabLimDetB(ami_mem, which, stldet);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetStabLimDetB");
}

void CVodeSolver::AMISetId(Model *model) { return; }

void CVodeSolver::AMISetSuppressAlg(bool flag) { return; }

void CVodeSolver::AMIReInit(realtype t0, N_Vector yy0, N_Vector yp0) {
    int status = CVodeReInit(ami_mem, t0, yy0);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeReInit");
}

void CVodeSolver::AMISensReInit(int ism, N_Vector *yS0, N_Vector *ypS0) {
    int status = CVodeSensReInit(ami_mem, ism, yS0);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSensReInit");
}

void CVodeSolver::AMISetSensParams(realtype *p, realtype *pbar, int *plist) {
    int status = CVodeSetSensParams(ami_mem, p, pbar, plist);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetSensParams");
}

void CVodeSolver::AMIGetDky(realtype t, int k, N_Vector dky) {
    int status = CVodeGetDky(ami_mem, t, k, dky);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeGetDky");
}

void CVodeSolver::AMIGetSens(realtype *tret, N_Vector *yySout) {
    int status = CVodeGetSens(ami_mem, tret, yySout);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeGetSens");
}
    
void CVodeSolver::AMIFree() {
    CVodeFree(&ami_mem);
    ami_mem = NULL;
}

void CVodeSolver::AMIAdjInit(long steps, int interp) {
    int status = CVodeAdjInit(ami_mem, steps, interp);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeAdjInit");
}

void CVodeSolver::AMICreateB(int lmm, int iter, int *which) {
    int status = CVodeCreateB(ami_mem, lmm, iter, which);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeCreateB");
}

void CVodeSolver::AMIReInitB(int which, realtype tB0, N_Vector yyB0,
                            N_Vector ypB0) {
    int status = CVodeReInitB(ami_mem, which, tB0, yyB0);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeReInitB");
}

void CVodeSolver::AMISStolerancesB(int which, realtype relTolB,
                                  realtype absTolB) {
    int status = CVodeSStolerancesB(ami_mem, which, relTolB, absTolB);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSStolerancesB");
}

void CVodeSolver::AMIQuadReInitB(int which, N_Vector yQB0) {
    int status = CVodeQuadReInitB(ami_mem, which, yQB0);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeQuadReInitB");
}

void CVodeSolver::AMIQuadSStolerancesB(int which, realtype reltolQB,
                                      realtype abstolQB) {
    int status = CVodeQuadSStolerancesB(ami_mem, which, reltolQB, abstolQB);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeQuadSStolerancesB");
}

int CVodeSolver::AMISolve(realtype tout, N_Vector yret, N_Vector ypret,
                          realtype *tret, int itask) {
    int status = CVode(ami_mem, tout, yret, tret, itask);
    if(status<0) {
        throw IntegrationFailure(status,*tret);
    } else{
        return status;
    }
}

int CVodeSolver::AMISolveF(realtype tout, N_Vector yret, N_Vector ypret,
                           realtype *tret, int itask, int *ncheckPtr) {
    int status = CVodeF(ami_mem, tout, yret, tret, itask, ncheckPtr);
    if(status<0) {
        throw IntegrationFailure(status,*tret);
    } else{
        return status;
    }
}

void CVodeSolver::AMISolveB(realtype tBout, int itaskB) {
    int status = CVodeB(ami_mem, tBout, itaskB);
    if(status != CV_SUCCESS)
         throw IntegrationFailure(status,tBout);
}

void CVodeSolver::AMISetMaxNumStepsB(int which, long mxstepsB) {
    int status = CVodeSetMaxNumStepsB(ami_mem, which, mxstepsB);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetMaxNumStepsB");
}

void CVodeSolver::AMIGetB(int which, realtype *tret, N_Vector yy, N_Vector yp) {
    int status = CVodeGetB(ami_mem, which, tret, yy);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeGetB");
}

void CVodeSolver::AMIGetQuadB(int which, realtype *tret, N_Vector qB) {
    int status = CVodeGetQuadB(ami_mem, which, tret, qB);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeGetQuadB");
}

void CVodeSolver::AMIDense(int nx) {
    int status = CVDense(ami_mem, nx);
    if(status != CV_SUCCESS)
        throw CvodeException(status,"CVDense");
}

void CVodeSolver::AMIDenseB(int which, int nx) {
    int status = CVDenseB(ami_mem, which, nx);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVDenseB");
}

void CVodeSolver::AMIBand(int nx, int ubw, int lbw) {
    int status = CVBand(ami_mem, nx, ubw, lbw);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVBand");
}

void CVodeSolver::AMIBandB(int which, int nx, int ubw, int lbw) {
    int status = CVBandB(ami_mem, which, nx, ubw, lbw);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVBandB");
}

void CVodeSolver::AMIDiag() {
    int status = CVDiag(ami_mem);
    if(status != CV_SUCCESS)
        throw CvodeException(status,"CVDiag");
}

void CVodeSolver::AMIDiagB(int which) {
    int status = CVDiagB(ami_mem, which);
        if(status != CV_SUCCESS)
            throw CvodeException(status,"CVDiagB");
}

void CVodeSolver::AMISpgmr(int prectype, int maxl) {
    int status = CVSpgmr(ami_mem, prectype, maxl);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVSpgmr");
}

void CVodeSolver::AMISpgmrB(int which, int prectype, int maxl) {
    int status = CVSpgmrB(ami_mem, which, prectype, maxl);
    if(status != CV_SUCCESS)
        throw CvodeException(status,"CVSpgmrB");
}

void CVodeSolver::AMISpbcg(int prectype, int maxl) {
    int status = CVSpbcg(ami_mem, prectype, maxl);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVSpbcg");
}

void CVodeSolver::AMISpbcgB(int which, int prectype, int maxl) {
    int status = CVSpbcgB(ami_mem, which, prectype, maxl);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVSpbcgB");
}

void CVodeSolver::AMISptfqmr(int prectype, int maxl) {
    int status = CVSptfqmr(ami_mem, prectype, maxl);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"AMISptfqmr");
}

void CVodeSolver::AMISptfqmrB(int which, int prectype, int maxl) {
    int status = CVSptfqmrB(ami_mem, which, prectype, maxl);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVSptfqmrB");
}

void CVodeSolver::AMIKLU(int nx, int nnz, int sparsetype) {
    int status = CVKLU(ami_mem, nx, nnz, sparsetype);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVKLU");
}

void CVodeSolver::AMIKLUSetOrdering(int ordering) {
    int status = CVKLUSetOrdering(ami_mem, ordering);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVKLUSetOrdering");
}

void CVodeSolver::AMIKLUSetOrderingB(int which, int ordering) {
    int status = CVKLUSetOrderingB(ami_mem, which, ordering);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVKLUSetOrderingB");
}

void CVodeSolver::AMIKLUB(int which, int nx, int nnz, int sparsetype) {
    int status = CVKLUB(ami_mem, which, nx, nnz, sparsetype);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVKLUB");
}

void CVodeSolver::AMIGetNumSteps(void *ami_mem, long *numsteps) {
    int status = CVodeGetNumSteps(ami_mem, numsteps);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeGetNumSteps");
}

void CVodeSolver::AMIGetNumRhsEvals(void *ami_mem, long *numrhsevals) {
    int status = CVodeGetNumRhsEvals(ami_mem, numrhsevals);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeGetNumRhsEvals");
}

void CVodeSolver::AMIGetNumErrTestFails(void *ami_mem, long *numerrtestfails) {
    int status = CVodeGetNumErrTestFails(ami_mem, numerrtestfails);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeGetNumErrTestFails");
}

void CVodeSolver::AMIGetNumNonlinSolvConvFails(void *ami_mem,
                                              long *numnonlinsolvconvfails) {
    int status = CVodeGetNumNonlinSolvConvFails(ami_mem, numnonlinsolvconvfails);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeGetNumNonlinSolvConvFails");
}

void CVodeSolver::AMIGetLastOrder(void *ami_ami_mem, int *order) {
    int status = CVodeGetLastOrder(ami_mem, order);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeGetLastOrder");
}

void *CVodeSolver::AMIGetAdjBmem(void *ami_mem, int which) {
    return CVodeGetAdjCVodeBmem(ami_mem, which);
}

void CVodeSolver::AMICalcIC(realtype tout1,TempData *tdata) { };

void CVodeSolver::AMICalcICB(int which, realtype tout1, N_Vector xB,
                             N_Vector dxB) {};

void CVodeSolver::AMISetStopTime(realtype tstop) {
    int status = CVodeSetStopTime(ami_mem, tstop);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetStopTime");
}

void CVodeSolver::turnOffRootFinding() {
    int status = CVodeRootInit(ami_mem, 0, NULL);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeRootInit");
}

CVodeSolver::~CVodeSolver() { AMIFree(); }

} // namespace amici
