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
#include <include/udata.h>

namespace amici {

void IDASolver::init(AmiVector *x, AmiVector *dx, realtype t) {
    int status = IDAInit(ami_mem, fxdot, RCONST(t), x->getNVector(), dx->getNVector());
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAInit");
}
void IDASolver::binit(int which, AmiVector *xB, AmiVector *dxB, realtype t) {
    int status = IDAInitB(ami_mem, which, fxBdot, RCONST(t), xB->getNVector(), dxB->getNVector());
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAInitB");
}
void IDASolver::qbinit(int which, AmiVector *qBdot) {
    int status = IDAQuadInitB(ami_mem, which, fqBdot, qBdot->getNVector());
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAQuadInitB");
}
void IDASolver::rootInit(int ne) {
    int status = IDARootInit(ami_mem, ne, froot);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDARootInit");
}
void IDASolver::sensInit1(AmiVectorArray *sx, AmiVectorArray *sdx, const UserData *udata) {
    int status = IDASensInit(ami_mem, udata->nplist(), udata->sensmeth(), fsxdot,
                             sx->getNVectorArray(),sdx->getNVectorArray());
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASensInit");
}
void IDASolver::setDenseJacFn() {
    int status = IDADlsSetDenseJacFn(ami_mem, fJ);
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDADlsSetDenseJacFn");
}

void IDASolver::setSparseJacFn() {
    int status = IDASlsSetSparseJacFn(ami_mem, fJSparse);
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDASlsSetSparseJacFn");
}
void IDASolver::setBandJacFn() {
    int status = IDADlsSetBandJacFn(ami_mem, fJBand);
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDADlsSetBandJacFn");
}

void IDASolver::setJacTimesVecFn() {
    int status = IDASpilsSetJacTimesVecFn(ami_mem, fJv);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASpilsSetJacTimesVecFn");
}
void IDASolver::setDenseJacFnB(int which) {
    int status = IDADlsSetDenseJacFnB(ami_mem, which, fJB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDADlsSetDenseJacFnB");
}
void IDASolver::setSparseJacFnB(int which) {
    int status = IDASlsSetSparseJacFnB(ami_mem, which, fJSparseB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASlsSetSparseJacFnB");
}
void IDASolver::setBandJacFnB(int which) {
    int status = IDADlsSetBandJacFnB(ami_mem, which, fJBandB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDADlsSetBandJacFnB");
}
void IDASolver::setJacTimesVecFnB(int which) {
    int status = IDASpilsSetJacTimesVecFnB(ami_mem, which, fJvB);
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

    AmiVector id(model->idlist);
    
    int status = IDASetId(ami_mem, id.getNVector());
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDASetMaxNumSteps");

}

void IDASolver::AMISetSuppressAlg(bool flag) {
    int status = IDASetSuppressAlg(ami_mem, flag);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetSuppressAlg");
}
void IDASolver::AMIReInit(realtype t0, AmiVector *yy0, AmiVector *yp0) {
    int status = IDAReInit(ami_mem, t0, yy0->getNVector(), yp0->getNVector());
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAReInit");
}
void IDASolver::AMISensReInit(int ism, AmiVectorArray *yS0, AmiVectorArray *ypS0) {
    int status = IDASensReInit(ami_mem, ism, yS0->getNVectorArray(), ypS0->getNVectorArray());
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASensReInit");
}
void IDASolver::AMISetSensParams(realtype *p, realtype *pbar, int *plist) {
    int status = IDASetSensParams(ami_mem, p, pbar, plist);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetSensParams");
}
void IDASolver::AMIGetDky(realtype t, int k, AmiVector *dky) {
    int status = IDAGetDky(ami_mem, t, k, dky->getNVector());
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetDky");
}
void IDASolver::AMIGetSens(realtype *tret, AmiVectorArray *yySout) {
    int status = IDAGetSens(ami_mem, tret, yySout->getNVectorArray());
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
void IDASolver::AMIReInitB(int which, realtype tB0, AmiVector *yyB0,
                          AmiVector *ypB0) {
    int status = IDAReInitB(ami_mem, which, tB0, yyB0->getNVector(), ypB0->getNVector());
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAReInitB");
}
void IDASolver::AMISStolerancesB(int which, realtype relTolB, realtype absTolB) {
    int status = IDASStolerancesB(ami_mem, which, relTolB, absTolB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASStolerancesB");
}
void IDASolver::AMIQuadReInitB(int which, AmiVector *yQB0) {
    int status = IDAQuadReInitB(ami_mem, which, yQB0->getNVector());
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAQuadReInitB");
}
void IDASolver::AMIQuadSStolerancesB(int which, realtype reltolQB,
                                    realtype abstolQB) {
    int status = IDAQuadSStolerancesB(ami_mem, which, reltolQB, abstolQB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAQuadSStolerancesB");
}
int IDASolver::AMISolve(realtype tout, AmiVector *yret, AmiVector *ypret,
                        realtype *tret, int itask) {
    int status = IDASolve(ami_mem, tout, tret, yret->getNVector(), ypret->getNVector(), itask);
    if(status<0) {
        throw IntegrationFailure(status,*tret);
    } else{
        return status;
    }
}
int IDASolver::AMISolveF(realtype tout, AmiVector *yret, AmiVector *ypret,
                         realtype *tret, int itask, int *ncheckPtr) {
    int status = IDASolveF(ami_mem, tout, tret, yret->getNVector(), ypret->getNVector(), itask, ncheckPtr);
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
void IDASolver::AMIGetB(int which, realtype *tret, AmiVector *yy, AmiVector *yp) {
    int status = IDAGetB(ami_mem, which, tret, yy->getNVector(), yp->getNVector());
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetB");
}
void IDASolver::AMIGetQuadB(int which, realtype *tret, AmiVector *qB) {
    int status = IDAGetQuadB(ami_mem, which, tret, qB->getNVector());
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

void IDASolver::AMICalcIC(realtype tout1, AmiVector *x, AmiVector *dx) {
    int status = IDACalcIC(ami_mem, IDA_YA_YDP_INIT, tout1);
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDACalcIC");
    status = IDAGetConsistentIC(ami_mem, x->getNVector(), dx->getNVector());
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDACalcIC");
}

void IDASolver::AMICalcICB(int which, realtype tout1, AmiVector *xB,
                          AmiVector *dxB) {
    int status = IDACalcICB(ami_mem, which, tout1, xB->getNVector(), dxB->getNVector());
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
    
    /** Jacobian of xdot with respect to states x
     * @param[in] N number of state variables
     * @param[in] t timepoint
     * @param[in] cj scaling factor, inverse of the step size
     * @param[in] x Vector with the states
     * @param[in] dx Vector with the derivative states
     * @param[in] xdot Vector with the right hand side
     * @param[out] J Matrix to which the Jacobian will be written
     * @param[in] user_data object with user input @type UserData
     * @param[in] tmp1 temporary storage vector
     * @param[in] tmp2 temporary storage vector
     * @param[in] tmp3 temporary storage vector
     * @return status flag indicating successful execution
     **/
    int IDASolver::fJ(long int N, realtype t, realtype cj, N_Vector x, N_Vector dx,
                  N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1,
                  N_Vector tmp2, N_Vector tmp3) {
        Model_DAE *model = (Model_DAE*) user_data;
        model->fdwdx(t,x);
        memset(J->data,0.0,sizeof(realtype)*N);
        model->model_J(J->data,t,N_VGetArrayPointer(x),model->p.data(),model->k.data(),model->h.data(),
                cj,N_VGetArrayPointer(dx),model->w.data(),model->dwdx.data());
        return checkVals(N,J->data,"Jacobian");
    }
    
    /** Jacobian of xBdot with respect to adjoint state xB
     * @param[in] NeqBdot number of adjoint state variables
     * @param[in] t timepoint
     * @param[in] cj scaling factor, inverse of the step size
     * @param[in] x Vector with the states
     * @param[in] dx Vector with the derivative states
     * @param[in] xB Vector with the adjoint states
     * @param[in] dxB Vector with the adjoint derivative states
     * @param[in] xBdot Vector with the adjoint right hand side
     * @param[out] JB Matrix to which the Jacobian will be written
     * @param[in] user_data object with user input @type UserData
     * @param[in] tmp1B temporary storage vector
     * @param[in] tmp2B temporary storage vector
     * @param[in] tmp3B temporary storage vector
     * @return status flag indicating successful execution
     **/
    int IDASolver::fJB(long int NeqBdot, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB,
                   N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B,
                   N_Vector tmp2B, N_Vector tmp3B) {
        Model_DAE *model = (Model_DAE*) user_data;
        model->fdwdx(t,x);
        memset(JB->data,0.0,sizeof(realtype)*NeqBdot);
        if(model->model_JB(JB->data,t,N_VGetArrayPointer(x),model->p.data(),model->k.data(),model->h.data(),
                     cj,N_VGetArrayPointer(xB),N_VGetArrayPointer(dx),N_VGetArrayPointer(dxB),
                     model->w.data(),model->dwdx.data()) != AMICI_SUCCESS)
        return AMICI_ERROR;
        return checkVals(NeqBdot,JB->data,"Jacobian");
    }
    
    /** J in sparse form (for sparse solvers from the SuiteSparse Package)
     * @param[in] t timepoint
     * @param[in] cj scalar in Jacobian (inverse stepsize)
     * @param[in] x Vector with the states
     * @param[in] dx Vector with the derivative states
     *N_Vector
     * @param[in] xdot Vector with the right hand side
     * @param[out] J Matrix to which the Jacobian will be written
     * @param[in] user_data object with user input @type UserData
     * @param[in] tmp1 temporary storage vector
     * @param[in] tmp2 temporary storage vector
     * @param[in] tmp3 temporary storage vector
     * @return status flag indicating successful execution
     */
    int IDASolver::fJSparse(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, SlsMat J,
                        void *user_data, N_Vector tmp1, N_Vector tmp2,
                        N_Vector tmp3) {
        Model_DAE *model = (Model_DAE*) user_data;
        model->fdwdx(t,x);
        memset(J->data,0.0,sizeof(realtype)*J->NNZ);
        model->model_JSparse(J->data,t,N_VGetArrayPointer(x),model->p.data(),model->k.data(),model->h.data(),
                      cj,N_VGetArrayPointer(dx),model->w.data(),model->dwdx.data());
        return checkVals(J->NNZ,J->data,"Jacobian");
    }

    /** JB in sparse form (for sparse solvers from the SuiteSparse Package)
     * @param[in] t timepoint
     * @param[in] cj scalar in Jacobian
     * @param[in] x Vector with the states
     * @param[in] dx Vector with the derivative states
     * @param[in] xB Vector with the adjoint states
     * @param[in] dxB Vector with the adjoint derivative states
     * @param[in] xBdot Vector with the adjoint right hand side
     * @param[out] JB Matrix to which the Jacobian will be written
     * @param[in] user_data object with user input @type UserData
     * @param[in] tmp1B temporary storage vector
     * @param[in] tmp2B temporary storage vector
     * @param[in] tmp3B temporary storage vector
     * @return status flag indicating successful execution
     */
    int IDASolver::fJSparseB(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot,
                         SlsMat JB, void *user_data, N_Vector tmp1B,
                         N_Vector tmp2B, N_Vector tmp3B) {
        Model_DAE *model = (Model_DAE*) user_data;
        model->fdwdx(t,x);
        memset(JB->data,0.0,sizeof(realtype)*JB->NNZ);
        if(model->model_JSparseB(JB->data,t,N_VGetArrayPointer(x),model->p.data(),model->k.data(),model->h.data(),
                           cj,N_VGetArrayPointer(xB),N_VGetArrayPointer(dx),N_VGetArrayPointer(dxB),
                           model->w.data(),model->dwdx.data()) != AMICI_SUCCESS)
            return AMICI_ERROR;
        return checkVals(JB->NNZ,JB->data,"Jacobian");
    }
    
    /** J in banded form (for banded solvers)
     * @param[in] N number of states
     * @param[in] mupper upper matrix bandwidth
     * @param[in] mlower lower matrix bandwidth
     * @param[in] t timepoint
     * @param[in] cj scalar in Jacobian (inverse stepsize)
     * @param[in] x Vector with the states
     * @param[in] dx Vector with the derivative states
     * @param[in] xdot Vector with the right hand side
     * @param[out] J Matrix to which the Jacobian will be written
     * @param[in] user_data object with user input @type UserData
     * @param[in] tmp1 temporary storage vector
     * @param[in] tmp2 temporary storage vector
     * @param[in] tmp3 temporary storage vector
     * @return status flag indicating successful execution
     */
    int IDASolver::fJBand(long int N, long int mupper, long int mlower, realtype t, realtype cj,
                      N_Vector x, N_Vector dx, N_Vector xdot, DlsMat J, void *user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
        return fJ(N,t,cj,x,dx,xdot,J,user_data,tmp1,tmp2,tmp3);
    }
    
    /** JB in banded form (for banded solvers)
     * @param[in] NeqBdot number of states
     * @param[in] mupper upper matrix bandwidth
     * @param[in] mlower lower matrix bandwidth
     * @param[in] t timepoint
     * @param[in] cj scalar in Jacobian (inverse stepsize)
     * @param[in] x Vector with the states
     * @param[in] dx Vector with the derivative states
     * @param[in] xB Vector with the adjoint states
     * @param[in] dxB Vector with the adjoint derivative states
     * @param[in] xBdot Vector with the adjoint right hand side
     * @param[out] JB Matrix to which the Jacobian will be written
     * @param[in] user_data object with user input @type UserData
     * @param[in] tmp1B temporary storage vector
     * @param[in] tmp2B temporary storage vector
     * @param[in] tmp3B temporary storage vector
     * @return status flag indicating successful execution
     */
    int IDASolver::fJBandB(long int NeqBdot, long int mupper, long int mlower,
                       realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot,
                       DlsMat JB, void *user_data, N_Vector tmp1B,
                       N_Vector tmp2B, N_Vector tmp3B) {
        return fJB(NeqBdot,t,cj,x,dx,xB,dxB,xBdot,JB,user_data,tmp1B,tmp2B,tmp3B);
    }
    

    
    /** Matrix vector product of J with a vector v (for iterative solvers)
     * @param[in] t timepoint @type realtype
     * @param[in] cj scaling factor, inverse of the step size
     * @param[in] x Vector with the states
     * @param[in] dx Vector with the derivative states
     * @param[in] xdot Vector with the right hand side
     * @param[in] v Vector with which the Jacobian is multiplied
     * @param[out] Jv Vector to which the Jacobian vector product will be
     *written
     * @param[in] user_data object with user input @type UserData
     * @param[in] tmp1 temporary storage vector
     * @param[in] tmp2 temporary storage vector
     * @return status flag indicating successful execution
     **/
    int IDASolver::fJv(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, N_Vector v, N_Vector Jv,
                   realtype cj, void *user_data, N_Vector tmp1, N_Vector tmp2) {
        Model_DAE *model = (Model_DAE*) user_data;
        model->fdwdx(t,x);
        memset(N_VGetArrayPointer(Jv),0.0,sizeof(realtype)*model->nx);
        if(model->model_Jv(N_VGetArrayPointer(Jv),t,N_VGetArrayPointer(x),model->p.data(),model->k.data(),model->h.data(),
                     cj,N_VGetArrayPointer(dx),N_VGetArrayPointer(v),model->w.data(),model->dwdx.data()) != AMICI_SUCCESS)
            return AMICI_ERROR;
        return checkVals(model->nx,N_VGetArrayPointer(Jv),"Jacobian");
    }
    
    /** Matrix vector product of JB with a vector v (for iterative solvers)
     * @param[in] t timepoint
     * @param[in] x Vector with the states
     * @param[in] dx Vector with the derivative states
     * @param[in] xB Vector with the adjoint states
     * @param[in] dxB Vector with the adjoint derivative states
     * @type N_Vector
     * @param[in] xBdot Vector with the adjoint right hand side
     * @param[in] vB Vector with which the Jacobian is multiplied
     * @param[out] JvB Vector to which the Jacobian vector product will be
     *written
     * @param[in] cj scalar in Jacobian (inverse stepsize)
     * @param[in] user_data object with user input @type UserData
     * @param[in] tmpB1 temporary storage vector
     * @param[in] tmpB2 temporary storage vector
     * @return status flag indicating successful execution
     **/
    int IDASolver::fJvB(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot,
                    N_Vector vB, N_Vector JvB, realtype cj, void *user_data,
                    N_Vector tmpB1, N_Vector tmpB2) {
        Model_DAE *model = (Model_DAE*) user_data;
        model->fdwdx(t,x);
        memset(N_VGetArrayPointer(JvB),0.0,sizeof(realtype)*model->nx);
        if(model->model_JvB(N_VGetArrayPointer(JvB),t,N_VGetArrayPointer(x),model->p.data(),model->k.data(),model->h.data(),
                      cj,N_VGetArrayPointer(xB),N_VGetArrayPointer(dx),N_VGetArrayPointer(dxB),
                      N_VGetArrayPointer(vB),model->w.data(),model->dwdx.data()) != AMICI_SUCCESS)
            return AMICI_ERROR;
        return checkVals(model->nx,N_VGetArrayPointer(JvB),"Jacobian");
    }
    
    /** Event trigger function for events
     * @param[in] t timepoint
     * @param[in] x Vector with the states
     * @param[in] dx Vector with the derivative states
     * @param[out] root array with root function values
     * @param[in] user_data object with user input @type UserData
     * @return status flag indicating successful execution
     */
    int IDASolver::froot(realtype t, N_Vector x, N_Vector dx, realtype *root,
                     void *user_data) {
        Model_DAE *model = (Model_DAE*) user_data;
        memset(root,0.0,sizeof(realtype)*model->ne);
        if(model->model_root(root,t,N_VGetArrayPointer(x),model->p.data(),model->k.data(),model->h.data(),
                       N_VGetArrayPointer(dx)) != AMICI_SUCCESS)
            return AMICI_ERROR;
        return checkVals(model->ne,root,"root function");
    }
    
    /** residual function of the DAE
     * @param[in] t timepoint
     * @param[in] x Vector with the states
     * @param[in] dx Vector with the derivative states
     * @param[out] xdot Vector with the right hand side
     * @param[in] user_data object with user input @type UserData
     * @return status flag indicating successful execution
     */
    int IDASolver::fxdot(realtype t, N_Vector x, N_Vector dx, N_Vector xdot,
                     void *user_data) {
        Model_DAE *model = (Model_DAE*) user_data;
        model->fw(t,x);
        memset(N_VGetArrayPointer(xdot),0.0,sizeof(realtype)*model->nx);
        model->model_xdot(N_VGetArrayPointer(xdot),t,N_VGetArrayPointer(x),model->p.data(),model->k.data(),model->h.data(),
                   N_VGetArrayPointer(dx),model->w.data());
        return checkVals(model->nx,N_VGetArrayPointer(xdot),"residual function");
    }
    
    /** Right hand side of differential equation for adjoint state xB
     * @param[in] t timepoint
     * @param[in] x Vector with the states
     * @param[in] dx Vector with the derivative states
     * @param[in] xB Vector with the adjoint states
     * @param[in] dxB Vector with the adjoint derivative states
     * @param[out] xBdot Vector with the adjoint right hand side
     * @param[in] user_data object with user input @type UserData
     * @return status flag indicating successful execution
     */
    int IDASolver::fxBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB,
                      N_Vector dxB, N_Vector xBdot, void *user_data) {
        Model_DAE *model = (Model_DAE*) user_data;
        model->fdwdx(t,x);
        memset(N_VGetArrayPointer(xBdot),0.0,sizeof(realtype)*model->nx);
        model->model_xBdot(N_VGetArrayPointer(xBdot),t,N_VGetArrayPointer(x),model->p.data(),model->k.data(),model->h.data(),
                    N_VGetArrayPointer(xB),N_VGetArrayPointer(dx),N_VGetArrayPointer(dxB),
                    model->w.data(),model->dwdx.data());
        return checkVals(model->nx,N_VGetArrayPointer(xBdot),"adjoint residual function");
    }
    
    /** Right hand side of integral equation for quadrature states qB
     * @param[in] t timepoint @type realtype
     * @param[in] x Vector with the states @type N_Vector
     * @param[in] dx Vector with the derivative states (only DAE) @type
     *N_Vector
     * @param[in] xB Vector with the adjoint states @type N_Vector
     * @param[in] dxB Vector with the adjoint derivative states (only DAE)
     * @type N_Vector
     * @param[out] qBdot Vector with the adjoint quadrature right hand side
     * @type N_Vector
     * @param[in] user_data pointer to temp data object @type TempDat
     * @return status flag indicating successful execution @type int
     */
    int IDASolver::fqBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector qBdot,
                      void *user_data) {
        Model_DAE *model = (Model_DAE*) user_data;
        model->fdwdp(t,x);
        memset(N_VGetArrayPointer(qBdot),0.0,sizeof(realtype)*model->nJ*model->plist.size());
        realtype *qBdot_tmp = N_VGetArrayPointer(qBdot);
        for(int ip = 1; ip < model->plist.size(); ip++)
            if(model->model_qBdot(&qBdot_tmp[ip*model->nJ],model->plist[ip],t,N_VGetArrayPointer(x),model->p.data(),model->k.data(),model->h.data(),
                    N_VGetArrayPointer(xB),N_VGetArrayPointer(dx),N_VGetArrayPointer(dxB),
                    model->w.data(),model->dwdp.data()) != AMICI_SUCCESS)
                return AMICI_ERROR;
        return checkVals(model->nJ*model->plist.size(),N_VGetArrayPointer(qBdot),"adjoint quadrature function");
    }
    
    /** Right hand side of differential equation for state sensitivities sx
     * @param[in] Ns number of parameters
     * @param[in] t timepoint
     * @param[in] x Vector with the states
     * @param[in] dx Vector with the derivative states
     * @param[in] xdot Vector with the right hand side
     * @param[in] ip parameter index
     * @param[in] sx Vector with the state sensitivities
     * @param[in] sdx Vector with the derivative state sensitivities
     * @param[out] sxdot Vector with the sensitivity right hand side
     * @param[in] user_data object with user input @type UserData
     * @param[in] tmp1 temporary storage vector
     * @param[in] tmp2 temporary storage vector
     * @param[in] tmp3 temporary storage vector
     * @return status flag indicating successful execution
     */
    int IDASolver::fsxdot(int Ns, realtype t, N_Vector x, N_Vector dx, N_Vector xdot,
                      N_Vector *sx, N_Vector *sdx, N_Vector *sxdot, void *user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
        Model_DAE *model = (Model_DAE*) user_data;
        model->fM(t,x);
        model->fdxdotdp(t,x,dx);
        fJSparse(t,0.0,x,dx,nullptr,model->J,model,tmp1,tmp2,tmp3);// also calls dwdx & dx
        for(int ip = 0; ip < model->plist.size(); ip++){
            memset(N_VGetArrayPointer(sxdot[ip]),0.0,sizeof(realtype)*model->nx);
            if(model->model_sxdot(N_VGetArrayPointer(sxdot[ip]),t,N_VGetArrayPointer(x),model->p.data(),model->k.data(),model->h.data(),
                            model->plist[ip],N_VGetArrayPointer(dx),N_VGetArrayPointer(sx[ip]),N_VGetArrayPointer(sdx[ip]),
                            model->w.data(),model->dwdx.data(),model->M.data(),model->J->data,model->dxdotdp.data()) != AMICI_SUCCESS)
                return AMICI_ERROR;
            if(checkVals(model->nx,N_VGetArrayPointer(sxdot[ip]),"sensitivity rhs") != AMICI_SUCCESS)
                return AMICI_ERROR;
        }
        return AMICI_SUCCESS;
    }
                                            
IDASolver::~IDASolver() { AMIFree(); }

} // namespace amici
