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
    int status = IDAInit(ami_mem, model->fxdot, RCONST(t), x, dx);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAInit");
}
void IDASolver::binit(int which, N_Vector xB, N_Vector dxB, realtype t) {
    int status = IDAInitB(ami_mem, which, model->fxBdot, RCONST(t), xB, dxB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAInitB");
}
void IDASolver::qbinit(int which, N_Vector qBdot) {
    int status = IDAQuadInitB(ami_mem, which, fqBdot, qBdot);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAQuadInitB");
}
void IDASolver::rootInit(int ne) {
    int status = IDARootInit(ami_mem, ne, model->froot);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDARootInit");
}
void IDASolver::sensInit1(N_Vector *sx, N_Vector *sdx, const UserData *udata) {
    int status = IDASensInit(ami_mem, udata->nplist, udata->sensi_meth, fsxdot, sx,
                       sdx);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASensInit");
}
void IDASolver::setDenseJacFn() {
    int status = IDADlsSetDenseJacFn(ami_mem, model->fJ);
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDADlsSetDenseJacFn");
}

void IDASolver::setSparseJacFn() {
    int status = IDASlsSetSparseJacFn(ami_mem, model->fJSparse);
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDASlsSetSparseJacFn");
}
void IDASolver::setBandJacFn() {
    int status = IDADlsSetBandJacFn(ami_mem, model->fJBand);
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDADlsSetBandJacFn");
}

void IDASolver::setJacTimesVecFn() {
    int status = IDASpilsSetJacTimesVecFn(ami_mem, model->fJv);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASpilsSetJacTimesVecFn");
}
void IDASolver::setDenseJacFnB(int which) {
    int status = IDADlsSetDenseJacFnB(ami_mem, which, model->fJB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDADlsSetDenseJacFnB");
}
void IDASolver::setSparseJacFnB(int which) {
    int status = IDASlsSetSparseJacFnB(ami_mem, which, model->fJSparseB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASlsSetSparseJacFnB");
}
void IDASolver::setBandJacFnB(int which) {
    int status = IDADlsSetBandJacFnB(ami_mem, which, model->fJBandB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDADlsSetBandJacFnB");
}
void IDASolver::setJacTimesVecFnB(int which) {
    int status = IDASpilsSetJacTimesVecFnB(ami_mem, which, model->fJvB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASpilsSetJacTimesVecFnB");
}
void *IDASolver::AMICreate(int lmm, int iter) {
    return IDACreate();
}

void IDASolver::AMISStolerances(double rtol, double atol) {
    int status = IDASStolerances(ami_mem, rtol, atol);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASStolerances");
}
void IDASolver::AMISensEEtolerances() {
    int status = IDASensEEtolerances(ami_mem);
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDASensEEtolerances");
}

void IDASolver::AMISetSensErrCon(bool error_corr) {
    int status = IDASetSensErrCon(ami_mem, error_corr);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetSensErrCon");
}
void IDASolver::AMISetQuadErrConB(int which, bool flag) {
    int status = IDASetQuadErrConB(ami_mem, which, flag);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetQuadErrConB");
}
void IDASolver::AMIGetRootInfo(int *rootsfound) {
    int status = IDAGetRootInfo(ami_mem, rootsfound);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetRootInfo");
}
void IDASolver::AMISetErrHandlerFn() {
    int status = IDASetErrHandlerFn(ami_mem, wrapErrHandlerFn, NULL);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetErrHandlerFn");
}
void IDASolver::AMISetUserData(void *user_data) {
    int status = IDASetUserData(ami_mem, user_data);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetUserData");
}
void IDASolver::AMISetUserDataB(int which, void *user_data) {
    int status = IDASetUserDataB(ami_mem, which, user_data);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetUserDataB");
}
void IDASolver::AMISetMaxNumSteps(long mxsteps) {
    int status = IDASetMaxNumSteps(ami_mem, mxsteps);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetMaxNumSteps");
}
void IDASolver::AMISetStabLimDet(int stldet) {
}

void IDASolver::AMISetStabLimDetB(int which, int stldet) {
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

}

void IDASolver::AMISetSuppressAlg(bool flag) {
    int status = IDASetSuppressAlg(ami_mem, flag);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetSuppressAlg");
}
void IDASolver::AMIReInit(realtype t0, N_Vector yy0, N_Vector yp0) {
    int status = IDAReInit(ami_mem, t0, yy0, yp0);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAReInit");
}
void IDASolver::AMISensReInit(int ism, N_Vector *yS0, N_Vector *ypS0) {
    int status = IDASensReInit(ami_mem, ism, yS0, ypS0);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASensReInit");
}
void IDASolver::AMISetSensParams(realtype *p, realtype *pbar, int *plist) {
    int status = IDASetSensParams(ami_mem, p, pbar, plist);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetSensParams");
}
void IDASolver::AMIGetDky(realtype t, int k, N_Vector dky) {
    int status = IDAGetDky(ami_mem, t, k, dky);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetDky");
}
void IDASolver::AMIGetSens(realtype *tret, N_Vector *yySout) {
    int status = IDAGetSens(ami_mem, tret, yySout);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetSens");
}
void IDASolver::AMIFree() {
    IDAFree(&ami_mem);
    ami_mem = NULL;
}

void IDASolver::AMIAdjInit(long steps, int interp) {
    int status = IDAAdjInit(ami_mem, steps, interp);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAAdjInit");
}
void IDASolver::AMICreateB(int lmm, int iter, int *which) {
    int status = IDACreateB(ami_mem, which);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDACreateB");
}
void IDASolver::AMIReInitB(int which, realtype tB0, N_Vector yyB0,
                          N_Vector ypB0) {
    int status = IDAReInitB(ami_mem, which, tB0, yyB0, ypB0);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAReInitB");
}
void IDASolver::AMISStolerancesB(int which, realtype relTolB, realtype absTolB) {
    int status = IDASStolerancesB(ami_mem, which, relTolB, absTolB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASStolerancesB");
}
void IDASolver::AMIQuadReInitB(int which, N_Vector yQB0) {
    int status = IDAQuadReInitB(ami_mem, which, yQB0);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAQuadReInitB");
}
void IDASolver::AMIQuadSStolerancesB(int which, realtype reltolQB,
                                    realtype abstolQB) {
    int status = IDAQuadSStolerancesB(ami_mem, which, reltolQB, abstolQB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAQuadSStolerancesB");
}
int IDASolver::AMISolve(realtype tout, N_Vector yret, N_Vector ypret,
                        realtype *tret, int itask) {
    int status = IDASolve(ami_mem, tout, tret, yret, ypret, itask);
    if(status<0) {
        throw IntegrationFailure(status,*tret);
    } else{
        return status;
    }
}
int IDASolver::AMISolveF(realtype tout, N_Vector yret, N_Vector ypret,
                         realtype *tret, int itask, int *ncheckPtr) {
    int status = IDASolveF(ami_mem, tout, tret, yret, ypret, itask, ncheckPtr);
    if(status<0) {
        throw IntegrationFailure(status,*tret);
    } else{
        return status;
    }
}
void IDASolver::AMISolveB(realtype tBout, int itaskB) {
    int status = IDASolveB(ami_mem, tBout, itaskB);
    if(status != IDA_SUCCESS)
         throw IntegrationFailure(status,tBout);
}
void IDASolver::AMISetMaxNumStepsB(int which, long mxstepsB) {
    int status = IDASetMaxNumStepsB(ami_mem, which, mxstepsB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetMaxNumStepsB");
}
void IDASolver::AMIGetB(int which, realtype *tret, N_Vector yy, N_Vector yp) {
    int status = IDAGetB(ami_mem, which, tret, yy, yp);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetB");
}
void IDASolver::AMIGetQuadB(int which, realtype *tret, N_Vector qB) {
    int status = IDAGetQuadB(ami_mem, which, tret, qB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetQuadB");
}
void IDASolver::AMIDense(int nx) {
    int status = IDADense(ami_mem, nx);
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDADense");
}

void IDASolver::AMIDenseB(int which, int nx) {
    int status = IDADenseB(ami_mem, which, nx);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDADenseB");
}
                                            
void IDASolver::AMIBand(int nx, int ubw, int lbw) {
    int status = IDABand(ami_mem, nx, ubw, lbw);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDABand");
}
void IDASolver::AMIBandB(int which, int nx, int ubw, int lbw) {
    int status = IDABandB(ami_mem, which, nx, ubw, lbw);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDABandB");
}
void IDASolver::AMIDiag() {
    throw AmiException("Diag Solver was not implemented for DAEs");
}

void IDASolver::AMIDiagB(int which) {
    throw AmiException("Diag Solver was not implemented for DAEs");
}

void IDASolver::AMISpgmr(int prectype, int maxl) {
    int status = IDASpgmr(ami_mem, maxl);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASpgmr");
}
void IDASolver::AMISpgmrB(int which, int prectype, int maxl) {
    int status = IDASpgmrB(ami_mem, which, maxl);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASpgmrB");
}
void IDASolver::AMISpbcg(int prectype, int maxl) {
    int status = IDASpbcg(ami_mem, maxl);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASpbcg");
}
void IDASolver::AMISpbcgB(int which, int prectype, int maxl) {
    int status = IDASpbcgB(ami_mem, which, maxl);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASpbcgB");
}
void IDASolver::AMISptfqmr(int prectype, int maxl) {
    int status = IDASptfqmr(ami_mem, maxl);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASptfqmr");
}
void IDASolver::AMISptfqmrB(int which, int prectype, int maxl) {
    int status = IDASptfqmrB(ami_mem, which, maxl);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASptfqmrB");
}
void IDASolver::AMIKLU(int nx, int nnz, int sparsetype) {
    int status = IDAKLU(ami_mem, nx, nnz, sparsetype);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAKLU");
}
void IDASolver::AMIKLUSetOrdering(int ordering) {
    int status = IDAKLUSetOrdering(ami_mem, ordering);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAKLUSetOrdering");
}
void IDASolver::AMIKLUSetOrderingB(int which, int ordering) {
    int status = IDAKLUSetOrderingB(ami_mem, which, ordering);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAKLUSetOrderingB");
}
void IDASolver::AMIKLUB(int which, int nx, int nnz, int sparsetype) {
    int status = IDAKLUB(ami_mem, which, nx, nnz, sparsetype);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAKLUB");
}
void IDASolver::AMIGetNumSteps(void *ami_mem, long *numsteps) {
    int status = IDAGetNumSteps(ami_mem, numsteps);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetNumSteps");
}
void IDASolver::AMIGetNumRhsEvals(void *ami_mem, long *numrhsevals) {
    int status = IDAGetNumResEvals(ami_mem, numrhsevals);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetNumResEvals");
}
void IDASolver::AMIGetNumErrTestFails(void *ami_mem, long *numerrtestfails) {
    int status = IDAGetNumErrTestFails(ami_mem, numerrtestfails);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetNumErrTestFails");
}
void IDASolver::AMIGetNumNonlinSolvConvFails(void *ami_mem,
                                            long *numnonlinsolvconvfails) {
    int status = IDAGetNumNonlinSolvConvFails(ami_mem, numnonlinsolvconvfails);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetNumNonlinSolvConvFails");
}
void IDASolver::AMIGetLastOrder(void *ami_mem, int *order) {
    int status = IDAGetLastOrder(ami_mem, order);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetLastOrder");
}
void *IDASolver::AMIGetAdjBmem(void *ami_mem, int which) {
    return IDAGetAdjIDABmem(ami_mem, which);
}

void IDASolver::AMICalcIC(realtype tout1, TempData *tdata) {
    int status = IDACalcIC(ami_mem, IDA_YA_YDP_INIT, tout1);
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDACalcIC");
    status = IDAGetConsistentIC(ami_mem, tdata->x, tdata->dx);
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDACalcIC");
}

void IDASolver::AMICalcICB(int which, realtype tout1, N_Vector xB,
                          N_Vector dxB) {
    int status = IDACalcICB(ami_mem, which, tout1, xB, dxB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDACalcICB");
}
    
void IDASolver::AMISetStopTime(realtype tstop) {
    int status = IDASetStopTime(ami_mem, tstop);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetStopTime");
}
                                            
void IDASolver::turnOffRootFinding() {
    int status = IDARootInit(ami_mem, 0, NULL);
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDARootInit");
}

                                            
IDASolver::~IDASolver() { AMIFree(); }

} // namespace amici
