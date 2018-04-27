#include "amici/solver_idas.h"

#include "amici/misc.h"
#include "amici/model.h"
#include "amici/exception.h"
#include "amici/model_dae.h"

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
void IDASolver::sensInit1(AmiVectorArray *sx, AmiVectorArray *sdx, int nplist) {
    int status = IDASensInit(ami_mem, nplist, getSensitivityMethod(), fsxdot,
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
Solver *IDASolver::clone() const {
    return new IDASolver(*this);
}

void *IDASolver::AMICreate(int lmm, int iter) {
    return IDACreate();
}

void IDASolver::AMISStolerances(double rtol, double atol) {
    int status = IDASStolerances(ami_mem, rtol, atol);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASStolerances");
}
void IDASolver::AMISensSStolerances(double rtol, double *atol) {
    int status = IDASensSStolerances(ami_mem, rtol, atol);
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
void IDASolver::AMISetUserData(Model *model) {
    int status = IDASetUserData(ami_mem, model);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetUserData");
}
void IDASolver::AMISetUserDataB(int which, Model *model) {
    int status = IDASetUserDataB(ami_mem, which, model);
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

    N_Vector id = N_VMake_Serial(model->nx,const_cast<realtype*>(model->idlist.data()));
    
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
        solverWasCalled = true;
        return status;
    }
}
int IDASolver::AMISolveF(realtype tout, AmiVector *yret, AmiVector *ypret,
                         realtype *tret, int itask, int *ncheckPtr) {
    int status = IDASolveF(ami_mem, tout, tret, yret->getNVector(), ypret->getNVector(), itask, ncheckPtr);
    if(status<0) {
        throw IntegrationFailure(status,*tret);
    } else{
        solverWasCalled = true;
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
     * @param N number of state variables
     * @param t timepoint
     * @param cj scaling factor, inverse of the step size
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xdot Vector with the right hand side
     * @param J Matrix to which the Jacobian will be written
     * @param user_data object with user input @type Model_DAE
     * @param tmp1 temporary storage vector
     * @param tmp2 temporary storage vector
     * @param tmp3 temporary storage vector
     * @return status flag indicating successful execution
     **/
    int IDASolver::fJ(long int N, realtype t, realtype cj, N_Vector x, N_Vector dx,
                  N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1,
                  N_Vector tmp2, N_Vector tmp3) {
        Model_DAE *model = static_cast<Model_DAE*>(user_data);
        model->fJ(t,cj, x, dx, xdot, J);
        return model->checkFinite(N,J->data,"Jacobian");
    }
    
    /** Jacobian of xBdot with respect to adjoint state xB
     * @param NeqBdot number of adjoint state variables
     * @param t timepoint
     * @param cj scaling factor, inverse of the step size
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xB Vector with the adjoint states
     * @param dxB Vector with the adjoint derivative states
     * @param xBdot Vector with the adjoint right hand side
     * @param JB Matrix to which the Jacobian will be written
     * @param user_data object with user input @type Model_DAE
     * @param tmp1B temporary storage vector
     * @param tmp2B temporary storage vector
     * @param tmp3B temporary storage vector
     * @return status flag indicating successful execution
     **/
    int IDASolver::fJB(long int NeqBdot, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB,
                   N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B,
                   N_Vector tmp2B, N_Vector tmp3B) {
        Model_DAE *model = static_cast<Model_DAE*>(user_data);
        model->fJB(t, cj, x, dx, xB, dxB, JB);
        return model->checkFinite(NeqBdot,JB->data,"Jacobian");
    }
    
    /** J in sparse form (for sparse solvers from the SuiteSparse Package)
     * @param t timepoint
     * @param cj scalar in Jacobian (inverse stepsize)
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xdot Vector with the right hand side
     * @param J Matrix to which the Jacobian will be written
     * @param user_data object with user input @type Model_DAE
     * @param tmp1 temporary storage vector
     * @param tmp2 temporary storage vector
     * @param tmp3 temporary storage vector
     * @return status flag indicating successful execution
     */
    int IDASolver::fJSparse(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, SlsMat J,
                        void *user_data, N_Vector tmp1, N_Vector tmp2,
                        N_Vector tmp3) {
        Model_DAE *model = static_cast<Model_DAE*>(user_data);
        model->fJSparse(t, cj, x, dx, J);
        return model->checkFinite(J->NNZ,J->data,"Jacobian");
    }

    /** JB in sparse form (for sparse solvers from the SuiteSparse Package)
     * @param t timepoint
     * @param cj scalar in Jacobian
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xB Vector with the adjoint states
     * @param dxB Vector with the adjoint derivative states
     * @param xBdot Vector with the adjoint right hand side
     * @param JB Matrix to which the Jacobian will be written
     * @param user_data object with user input @type Model_DAE
     * @param tmp1B temporary storage vector
     * @param tmp2B temporary storage vector
     * @param tmp3B temporary storage vector
     * @return status flag indicating successful execution
     */
    int IDASolver::fJSparseB(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot,
                         SlsMat JB, void *user_data, N_Vector tmp1B,
                         N_Vector tmp2B, N_Vector tmp3B) {
        Model_DAE *model = static_cast<Model_DAE*>(user_data);
        model->fJSparseB(t, cj, x, dx, xB, dxB, JB);
        return model->checkFinite(JB->NNZ,JB->data,"Jacobian");
    }
    
    /** J in banded form (for banded solvers)
     * @param N number of states
     * @param mupper upper matrix bandwidth
     * @param mlower lower matrix bandwidth
     * @param t timepoint
     * @param cj scalar in Jacobian (inverse stepsize)
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xdot Vector with the right hand side
     * @param J Matrix to which the Jacobian will be written
     * @param user_data object with user input @type Model_DAE
     * @param tmp1 temporary storage vector
     * @param tmp2 temporary storage vector
     * @param tmp3 temporary storage vector
     * @return status flag indicating successful execution
     */
    int IDASolver::fJBand(long int N, long int mupper, long int mlower, realtype t, realtype cj,
                      N_Vector x, N_Vector dx, N_Vector xdot, DlsMat J, void *user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
        return fJ(N,t,cj,x,dx,xdot,J,user_data,tmp1,tmp2,tmp3);
    }
    
    /** JB in banded form (for banded solvers)
     * @param NeqBdot number of states
     * @param mupper upper matrix bandwidth
     * @param mlower lower matrix bandwidth
     * @param t timepoint
     * @param cj scalar in Jacobian (inverse stepsize)
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xB Vector with the adjoint states
     * @param dxB Vector with the adjoint derivative states
     * @param xBdot Vector with the adjoint right hand side
     * @param JB Matrix to which the Jacobian will be written
     * @param user_data object with user input @type Model_DAE
     * @param tmp1B temporary storage vector
     * @param tmp2B temporary storage vector
     * @param tmp3B temporary storage vector
     * @return status flag indicating successful execution
     */
    int IDASolver::fJBandB(long int NeqBdot, long int mupper, long int mlower,
                       realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot,
                       DlsMat JB, void *user_data, N_Vector tmp1B,
                       N_Vector tmp2B, N_Vector tmp3B) {
        return fJB(NeqBdot,t,cj,x,dx,xB,dxB,xBdot,JB,user_data,tmp1B,tmp2B,tmp3B);
    }
    

    
    /** Matrix vector product of J with a vector v (for iterative solvers)
     * @param t timepoint @type realtype
     * @param cj scaling factor, inverse of the step size
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xdot Vector with the right hand side
     * @param v Vector with which the Jacobian is multiplied
     * @param Jv Vector to which the Jacobian vector product will be
     *written
     * @param user_data object with user input @type Model_DAE
     * @param tmp1 temporary storage vector
     * @param tmp2 temporary storage vector
     * @return status flag indicating successful execution
     **/
    int IDASolver::fJv(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, N_Vector v, N_Vector Jv,
                   realtype cj, void *user_data, N_Vector tmp1, N_Vector tmp2) {
        Model_DAE *model = static_cast<Model_DAE*>(user_data);
        model->fJv(t, x, dx, v, Jv, cj);
        return model->checkFinite(model->nx,N_VGetArrayPointer(Jv),"Jacobian");
    }
    
    /** Matrix vector product of JB with a vector v (for iterative solvers)
     * @param t timepoint
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xB Vector with the adjoint states
     * @param dxB Vector with the adjoint derivative states
     * @param xBdot Vector with the adjoint right hand side
     * @param vB Vector with which the Jacobian is multiplied
     * @param JvB Vector to which the Jacobian vector product will be
     *written
     * @param cj scalar in Jacobian (inverse stepsize)
     * @param user_data object with user input @type Model_DAE
     * @param tmpB1 temporary storage vector
     * @param tmpB2 temporary storage vector
     * @return status flag indicating successful execution
     **/
    int IDASolver::fJvB(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot,
                    N_Vector vB, N_Vector JvB, realtype cj, void *user_data,
                    N_Vector tmpB1, N_Vector tmpB2) {
        Model_DAE *model = static_cast<Model_DAE*>(user_data);
        model->fJvB(t, x, dx, xB, dxB, vB, JvB, cj);
        return model->checkFinite(model->nx,N_VGetArrayPointer(JvB),"Jacobian");
    }
    
    /** Event trigger function for events
     * @param t timepoint
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param root array with root function values
     * @param user_data object with user input @type Model_DAE
     * @return status flag indicating successful execution
     */
    int IDASolver::froot(realtype t, N_Vector x, N_Vector dx, realtype *root,
                     void *user_data) {
        Model_DAE *model = static_cast<Model_DAE*>(user_data);
        model->froot(t,x,dx,root);
        return model->checkFinite(model->ne,root,"root function");
    }
    
    /** residual function of the DAE
     * @param t timepoint
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xdot Vector with the right hand side
     * @param user_data object with user input @type Model_DAE
     * @return status flag indicating successful execution
     */
    int IDASolver::fxdot(realtype t, N_Vector x, N_Vector dx, N_Vector xdot,
                     void *user_data) {
        Model_DAE *model = static_cast<Model_DAE*>(user_data);
        model->fxdot(t,x,dx,xdot);
        return model->checkFinite(model->nx,N_VGetArrayPointer(xdot),"residual function");
    }
    
    /** Right hand side of differential equation for adjoint state xB
     * @param t timepoint
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xB Vector with the adjoint states
     * @param dxB Vector with the adjoint derivative states
     * @param xBdot Vector with the adjoint right hand side
     * @param user_data object with user input @type Model_DAE
     * @return status flag indicating successful execution
     */
    int IDASolver::fxBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB,
                      N_Vector dxB, N_Vector xBdot, void *user_data) {
        Model_DAE *model = static_cast<Model_DAE*>(user_data);
        model->fxBdot(t, x, dx, xB, dxB, xBdot);
        return model->checkFinite(model->nx,N_VGetArrayPointer(xBdot),"adjoint residual function");
    }
    
    /** Right hand side of integral equation for quadrature states qB
     * @param t timepoint
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xB Vector with the adjoint states
     * @param dxB Vector with the adjoint derivative states
     * @param qBdot Vector with the adjoint quadrature right hand side
     * @param user_data pointer to temp data object @type Model_DAE
     * @return status flag indicating successful execution
     */
    int IDASolver::fqBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector qBdot,
                      void *user_data) {
        Model_DAE *model = static_cast<Model_DAE*>(user_data);
        model->fqBdot(t, x, dx, xB, dxB, qBdot);
        return model->checkFinite(model->nJ*model->nplist(),N_VGetArrayPointer(qBdot),"adjoint quadrature function");
    }
    
    /** Right hand side of differential equation for state sensitivities sx
     * @param Ns number of parameters
     * @param t timepoint
     * @param x Vector with the states
     * @param dx Vector with the derivative states
     * @param xdot Vector with the right hand side
     * @param sx Vector with the state sensitivities
     * @param sdx Vector with the derivative state sensitivities
     * @param sxdot Vector with the sensitivity right hand side
     * @param user_data object with user input @type Model_DAE
     * @param tmp1 temporary storage vector
     * @param tmp2 temporary storage vector
     * @param tmp3 temporary storage vector
     * @return status flag indicating successful execution
     */
    int IDASolver::fsxdot(int Ns, realtype t, N_Vector x, N_Vector dx, N_Vector xdot,
                      N_Vector *sx, N_Vector *sdx, N_Vector *sxdot, void *user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
        Model_DAE *model = static_cast<Model_DAE*>(user_data);
        for(int ip = 0; ip < model->nplist(); ip++){
            model->fsxdot(t, x, dx, ip, sx[ip], sdx[ip], sxdot[ip]);
            if(model->checkFinite(model->nx,N_VGetArrayPointer(sxdot[ip]),"sensitivity rhs") != AMICI_SUCCESS)
                return AMICI_RECOVERABLE_ERROR;
        }
        return AMICI_SUCCESS;
    }
                                            
IDASolver::~IDASolver() { AMIFree(); }

} // namespace amici
