#include "amici/misc.h"
#include "amici/model_dae.h"
#include "amici/solver_idas.h"
#include "amici/exception.h"

#include <idas/idas.h>
#include <idas/idas_impl.h>

#include "amici/sundials_linsol_wrapper.h"

#include <amd.h>
#include <btf.h>
#include <colamd.h>
#include <klu.h>

#define ONE     RCONST(1.0)

namespace amici {

void IDASolver::init(AmiVector *x, AmiVector *dx, realtype t) {
    int status = IDAInit(solverMemory.get(), fxdot, t, x->getNVector(), dx->getNVector());
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAInit");
}
void IDASolver::binit(int which, AmiVector *xB, AmiVector *dxB, realtype t) {
    int status = IDAInitB(solverMemory.get(), which, fxBdot, t, xB->getNVector(), dxB->getNVector());
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAInitB");
}
void IDASolver::qbinit(int which, AmiVector *qBdot) {
    int status = IDAQuadInitB(solverMemory.get(), which, fqBdot, qBdot->getNVector());
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAQuadInitB");
}
void IDASolver::rootInit(int ne) {
    int status = IDARootInit(solverMemory.get(), ne, froot);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDARootInit");
}
void IDASolver::sensInit1(AmiVectorArray *sx, AmiVectorArray *sdx, int nplist) {
    int status = IDASensInit(solverMemory.get(), nplist, static_cast<int>(getSensitivityMethod()), fsxdot,
                             sx->getNVectorArray(),sdx->getNVectorArray());
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASensInit");
}
void IDASolver::setDenseJacFn() {
    int status = IDASetJacFn(solverMemory.get(), fJ);
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDADlsSetDenseJacFn");
}

void IDASolver::setSparseJacFn() {
    int status = IDASetJacFn(solverMemory.get(), fJSparse);
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDASlsSetSparseJacFn");
}
void IDASolver::setBandJacFn() {
    int status = IDASetJacFn(solverMemory.get(), fJBand);
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDADlsSetBandJacFn");
}

void IDASolver::setJacTimesVecFn() {
    int status = IDASetJacTimes(solverMemory.get(), nullptr, fJv);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASpilsSetJacTimesVecFn");
}
void IDASolver::setDenseJacFnB(int which) {
    int status = IDASetJacFnB(solverMemory.get(), which, fJB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDADlsSetDenseJacFnB");
}
void IDASolver::setSparseJacFnB(int which) {
    int status = IDASetJacFnB(solverMemory.get(), which, fJSparseB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASlsSetSparseJacFnB");
}
void IDASolver::setBandJacFnB(int which) {
    int status = IDASetJacFnB(solverMemory.get(), which, fJBandB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDADlsSetBandJacFnB");
}
void IDASolver::setJacTimesVecFnB(int which) {
    int status = IDASetJacTimesB(solverMemory.get(), which, nullptr, fJvB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASpilsSetJacTimesVecFnB");
}
Solver *IDASolver::clone() const {
    return new IDASolver(*this);
}

void IDASolver::allocateSolver() {
    solverMemory = std::unique_ptr<void, std::function<void(void *)>>(IDACreate(),
                   [](void *ptr) { IDAFree(&ptr); });
}

void IDASolver::setSStolerances(double rtol, double atol) {
    int status = IDASStolerances(solverMemory.get(), rtol, atol);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASStolerances");
}
void IDASolver::setSensSStolerances(double rtol, double *atol) {
    int status = IDASensSStolerances(solverMemory.get(), rtol, atol);
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDASensEEtolerances");
}

void IDASolver::setSensErrCon(bool error_corr) {
    int status = IDASetSensErrCon(solverMemory.get(), error_corr);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetSensErrCon");
}
void IDASolver::setQuadErrConB(int which, bool flag) {
    int status = IDASetQuadErrConB(solverMemory.get(), which, flag);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetQuadErrConB");
}
void IDASolver::getRootInfo(int *rootsfound) const {
    int status = IDAGetRootInfo(solverMemory.get(), rootsfound);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetRootInfo");
}
void IDASolver::setErrHandlerFn() {
    int status = IDASetErrHandlerFn(solverMemory.get(), wrapErrHandlerFn, nullptr);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetErrHandlerFn");
}
void IDASolver::setUserData(Model *model) {
    int status = IDASetUserData(solverMemory.get(), model);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetUserData");
}
void IDASolver::setUserDataB(int which, Model *model) {
    int status = IDASetUserDataB(solverMemory.get(), which, model);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetUserDataB");
}
void IDASolver::setMaxNumSteps(long mxsteps) {
    int status = IDASetMaxNumSteps(solverMemory.get(), mxsteps);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetMaxNumSteps");
}
void IDASolver::setStabLimDet(int stldet) {
}

void IDASolver::setStabLimDetB(int which, int stldet) {
}

void IDASolver::setId(Model *model) {

    N_Vector id = N_VMake_Serial(model->nx_solver,const_cast<realtype*>(model->idlist.data()));

    int status = IDASetId(solverMemory.get(), id);
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDASetMaxNumSteps");

    N_VDestroy_Serial(id);
}

void IDASolver::setSuppressAlg(bool flag) {
    int status = IDASetSuppressAlg(solverMemory.get(), flag);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetSuppressAlg");
}
    
void IDASolver::resetState(void *ami_mem, N_Vector yy0, N_Vector yp0) {
    
    auto ida_mem = static_cast<IDAMem>(ami_mem);
    /* here we force the order in the next step to zero, and update the
     phi arrays, this is largely copied from IDAReInit with
     explanations from idas_impl.h
     */
    
    /* Initialize the phi array */
    
    N_VScale(ONE, yy0, ida_mem->ida_phi[0]);
    N_VScale(ONE, yp0, ida_mem->ida_phi[1]);
    
    /* Set step parameters */
    
    /* current order */
    ida_mem->ida_kk      = 0;
    
    /* Initial setup not done yet */
    
    ida_mem->ida_SetupDone = SUNFALSE;
}

void IDASolver::reInitPostProcessF(realtype *t, AmiVector *yout,
                                   AmiVector *ypout, realtype tnext) {
    reInitPostProcess(solverMemory.get(), t, yout, ypout,  tnext);
}

void IDASolver::reInitPostProcessB(int which, realtype *t, AmiVector *yBout,
                                   AmiVector *ypBout, realtype tnext) {
    reInitPostProcess(IDAGetAdjIDABmem(solverMemory.get(), which),
                      t, yBout, ypBout, tnext);
}

void IDASolver::reInitPostProcess(void *ami_mem, realtype *t,
                                  AmiVector *yout, AmiVector *ypout,
                                  realtype tout) {}

void IDASolver::reInit(realtype t0, AmiVector *yy0, AmiVector *yp0) {
    auto ida_mem = static_cast<IDAMem>(solverMemory.get());
    /* set time */
    ida_mem->ida_tn = t0;
    resetState(ida_mem, yy0->getNVector(), yp0->getNVector());
}
void IDASolver::sensReInit(AmiVectorArray *yS0, AmiVectorArray *ypS0) {
    auto ida_mem = static_cast<IDAMem>(solverMemory.get());
    /* Initialize znS[0] in the history array */
    
    for (int is=0; is<ida_mem->ida_Ns; is++)
        ida_mem->ida_cvals[is] = ONE;
    
    auto status = N_VScaleVectorArray(ida_mem->ida_Ns, ida_mem->ida_cvals,
                                     yS0->getNVectorArray(),
                                     ida_mem->ida_phiS[0]);
    if(status != IDA_SUCCESS)
        throw IDAException(IDA_VECTOROP_ERR,"IDASensReInit");
    
    status = N_VScaleVectorArray(ida_mem->ida_Ns, ida_mem->ida_cvals,
                                     ypS0->getNVectorArray(),
                                     ida_mem->ida_phiS[1]);
    if(status != IDA_SUCCESS)
        throw IDAException(IDA_VECTOROP_ERR,"IDASensReInit");
}
void IDASolver::setSensParams(realtype *p, realtype *pbar, int *plist) {
    int status = IDASetSensParams(solverMemory.get(), p, pbar, plist);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetSensParams");
}
void IDASolver::getDky(realtype t, int k, AmiVector *dky) const {
    int status = IDAGetDky(solverMemory.get(), t, k, dky->getNVector());
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetDky");
}
void IDASolver::getSens(realtype *tret, AmiVectorArray *yySout) const {
    int status = IDAGetSens(solverMemory.get(), tret, yySout->getNVectorArray());
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetSens");
}

void IDASolver::adjInit() {
    int status = IDAAdjInit(solverMemory.get(), static_cast<int>(maxsteps), static_cast<int>(interpType));
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAAdjInit");
}
void IDASolver::allocateSolverB(int *which) {
    int status = IDACreateB(solverMemory.get(), which);
    if (*which + 1 > static_cast<int>(solverMemoryB.size()))
        solverMemoryB.resize(*which + 1);
    solverMemoryB.at(*which) = std::unique_ptr<void, std::function<void(void *)>>
    (getAdjBmem(solverMemory.get(), *which), [](void *ptr){});
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDACreateB");
}
void IDASolver::reInitB(int which, realtype tB0, AmiVector *yyB0,
                          AmiVector *ypB0) {
    auto ida_memB = static_cast<IDAMem>(IDAGetAdjIDABmem(solverMemory.get(), which));
    ida_memB->ida_tn = tB0;
    resetState(ida_memB, yyB0->getNVector(), ypB0->getNVector());
}
void IDASolver::setSStolerancesB(int which, realtype relTolB, realtype absTolB) {
    int status = IDASStolerancesB(solverMemory.get(), which, relTolB, absTolB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASStolerancesB");
}
void IDASolver::quadReInitB(int which, AmiVector *yQB0) {
    auto ida_memB = static_cast<IDAMem>(IDAGetAdjIDABmem(solverMemory.get(), which));
    N_VScale(ONE, yQB0->getNVector(), ida_memB->ida_phiQ[0]);
}
void IDASolver::quadSStolerancesB(int which, realtype reltolQB,
                                    realtype abstolQB) {
    int status = IDAQuadSStolerancesB(solverMemory.get(), which, reltolQB, abstolQB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAQuadSStolerancesB");
}
int IDASolver::solve(realtype tout, AmiVector *yret, AmiVector *ypret,
                        realtype *tret, int itask) {
    int status = IDASolve(solverMemory.get(), tout, tret, yret->getNVector(), ypret->getNVector(), itask);
    if(status<0) {
        throw IntegrationFailure(status,*tret);
    }

    solverWasCalled = true;
    return status;
}
int IDASolver::solveF(realtype tout, AmiVector *yret, AmiVector *ypret,
                         realtype *tret, int itask, int *ncheckPtr) {
    int status = IDASolveF(solverMemory.get(), tout, tret, yret->getNVector(), ypret->getNVector(), itask, ncheckPtr);
    if(status<0) {
        throw IntegrationFailure(status,*tret);
    }

    solverWasCalled = true;
    return status;
}
void IDASolver::solveB(realtype tBout, int itaskB) {
    int status = IDASolveB(solverMemory.get(), tBout, itaskB);
    if(status != IDA_SUCCESS)
         throw IntegrationFailure(status,tBout);
}
void IDASolver::setMaxNumStepsB(int which, long mxstepsB) {
    int status = IDASetMaxNumStepsB(solverMemory.get(), which, mxstepsB);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetMaxNumStepsB");
}
void IDASolver::getB(int which, realtype *tret, AmiVector *yy, AmiVector *yp) const {
    int status = IDAGetB(solverMemory.get(), which, tret, yy->getNVector(), yp->getNVector());
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetB");
}
void IDASolver::getQuadB(int which, realtype *tret, AmiVector *qB) const {
    int status = IDAGetQuadB(solverMemory.get(), which, tret, qB->getNVector());
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetQuadB");
}

void IDASolver::diag() {
    throw AmiException("Diag Solver was not implemented for DAEs");
}

void IDASolver::diagB(int  /*which*/) {
    throw AmiException("Diag Solver was not implemented for DAEs");
}

void IDASolver::getNumSteps(void *ami_mem, long *numsteps) const {
    int status = IDAGetNumSteps(ami_mem, numsteps);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetNumSteps");
}
void IDASolver::getNumRhsEvals(void *ami_mem, long *numrhsevals) const {
    int status = IDAGetNumResEvals(ami_mem, numrhsevals);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetNumResEvals");
}
void IDASolver::getNumErrTestFails(void *ami_mem, long *numerrtestfails) const {
    int status = IDAGetNumErrTestFails(ami_mem, numerrtestfails);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetNumErrTestFails");
}
void IDASolver::getNumNonlinSolvConvFails(void *ami_mem,
                                            long *numnonlinsolvconvfails) const {
    int status = IDAGetNumNonlinSolvConvFails(ami_mem, numnonlinsolvconvfails);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetNumNonlinSolvConvFails");
}
void IDASolver::getLastOrder(void *ami_mem, int *order) const {
    int status = IDAGetLastOrder(ami_mem, order);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDAGetLastOrder");
}
void *IDASolver::getAdjBmem(void *ami_mem, int which) {
    return IDAGetAdjIDABmem(ami_mem, which);
}

void IDASolver::calcIC(realtype tout1, AmiVector *x, AmiVector *dx) {
    int status = IDACalcIC(solverMemory.get(), IDA_YA_YDP_INIT, tout1);
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDACalcIC");
    status = IDAGetConsistentIC(solverMemory.get(), x->getNVector(), dx->getNVector());
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDACalcIC");
}

void IDASolver::calcICB(int which, realtype tout1, AmiVector *xB,
                          AmiVector *dxB) {
    int status = IDACalcICB(solverMemory.get(), which, tout1, xB->getNVector(), dxB->getNVector());
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDACalcICB");
}

void IDASolver::setStopTime(realtype tstop) {
    int status = IDASetStopTime(solverMemory.get(), tstop);
    if(status != IDA_SUCCESS)
         throw IDAException(status,"IDASetStopTime");
}

void IDASolver::turnOffRootFinding() {
    int status = IDARootInit(solverMemory.get(), 0, nullptr);
    if(status != IDA_SUCCESS)
        throw IDAException(status,"IDARootInit");
}

int IDASolver::nplist() const {
    if (!solverMemory)
        throw AmiException("Solver has not been allocated, information is not available");
    auto IDA_mem = (IDAMem) solverMemory.get();
    return IDA_mem->ida_Ns;
}

int IDASolver::nx() const {
    if (!solverMemory)
        throw AmiException("Solver has not been allocated, information is not available");
    auto IDA_mem = (IDAMem) solverMemory.get();
    return NV_LENGTH_S(IDA_mem->ida_yy0);
}

const Model *IDASolver::getModel() const {
    if (!solverMemory)
        throw AmiException("Solver has not been allocated, information is not available");
    auto ida_mem = (IDAMem) solverMemory.get();
    return static_cast<Model *>(ida_mem->ida_user_data);
}

bool IDASolver::getMallocDone() const {
    if (!solverMemory)
        return false;
    auto ida_mem = (IDAMem) solverMemory.get();
    return ida_mem->ida_MallocDone;
}

bool IDASolver::getAdjMallocDone() const {
    if (!solverMemory)
        return false;
    auto ida_mem = (IDAMem) solverMemory.get();
    return ida_mem->ida_adjMallocDone;
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
int IDASolver::fJ(realtype t, realtype cj, N_Vector x, N_Vector dx,
                  N_Vector xdot, SUNMatrix J, void *user_data, N_Vector  /*tmp1*/,
                  N_Vector  /*tmp2*/, N_Vector  /*tmp3*/) {
    auto model = static_cast<Model_DAE *>(user_data);
    model->fJ(t, cj, x, dx, xdot, J);
    return model->checkFinite(SM_ROWS_D(J), SM_DATA_D(J), "Jacobian");
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
int IDASolver::fJB(realtype t, realtype cj, N_Vector x, N_Vector dx,
                   N_Vector xB, N_Vector dxB, N_Vector  /*xBdot*/, SUNMatrix JB,
                   void *user_data, N_Vector  /*tmp1B*/, N_Vector  /*tmp2B*/,
                   N_Vector  /*tmp3B*/) {
    auto model = static_cast<Model_DAE *>(user_data);
    model->fJB(t, cj, x, dx, xB, dxB, JB);
    return model->checkFinite(SM_ROWS_D(JB), SM_DATA_D(JB), "Jacobian");
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
int IDASolver::fJSparse(realtype t, realtype cj, N_Vector x, N_Vector dx,
                        N_Vector  /*xdot*/, SUNMatrix J, void *user_data,
                        N_Vector  /*tmp1*/, N_Vector  /*tmp2*/, N_Vector  /*tmp3*/) {
    auto model = static_cast<Model_DAE *>(user_data);
    model->fJSparse(t, cj, x, dx, J);
    return model->checkFinite(SM_NNZ_S(J), SM_DATA_S(J), "Jacobian");
}

void IDASolver::setLinearSolver()
{
    int status = IDASetLinearSolver(solverMemory.get(), linearSolver->get(), linearSolver->getMatrix());
    if(status != IDA_SUCCESS)
        throw IDAException(status, "setLinearSolver");

}

void IDASolver::setLinearSolverB(int which)
{
    int status = IDASetLinearSolverB(solverMemoryB[which].get(), which, linearSolverB->get(), linearSolverB->getMatrix());
    if(status != IDA_SUCCESS)
        throw IDAException(status, "setLinearSolverB");

}

void IDASolver::setNonLinearSolver()
{
    int status = IDASetNonlinearSolver(solverMemory.get(), nonLinearSolver->get());
    if(status != IDA_SUCCESS)
        throw CvodeException(status,"CVodeSetNonlinearSolver");
}

void IDASolver::setNonLinearSolverSens()
{
    if(getSensitivityOrder() < SensitivityOrder::first)
        return;
    if(getSensitivityMethod() != SensitivityMethod::forward)
        return;

    int status = IDA_SUCCESS;

    switch (ism) {
    case InternalSensitivityMethod::staggered:
        status = IDASetNonlinearSolverSensStg(solverMemory.get(), nonLinearSolverSens->get());
        break;
    case InternalSensitivityMethod::simultaneous:
        status = IDASetNonlinearSolverSensSim(solverMemory.get(), nonLinearSolverSens->get());
        break;
    case InternalSensitivityMethod::staggered1:
    default:
        throw AmiException("Unsupported internal sensitivity method selected: %d", ism);
    }

    if(status != IDA_SUCCESS)
        throw CvodeException(status,"CVodeSolver::setNonLinearSolverSens");
}

void IDASolver::setNonLinearSolverB(int which)
{
    int status = IDASetNonlinearSolverB(solverMemory.get(), which, nonLinearSolverB->get());
    if(status != IDA_SUCCESS)
        throw CvodeException(status,"CVodeSetNonlinearSolverB");
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
int IDASolver::fJSparseB(realtype t, realtype cj, N_Vector x, N_Vector dx,
                         N_Vector xB, N_Vector dxB, N_Vector  /*xBdot*/,
                         SUNMatrix JB, void *user_data, N_Vector  /*tmp1B*/,
                         N_Vector  /*tmp2B*/, N_Vector  /*tmp3B*/) {
    auto model = static_cast<Model_DAE *>(user_data);
    model->fJSparseB(t, cj, x, dx, xB, dxB, JB);
    return model->checkFinite(SM_NNZ_S(JB), SM_DATA_S(JB), "Jacobian");
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
int IDASolver::fJBand(realtype t, realtype cj, N_Vector x, N_Vector dx,
                      N_Vector xdot, SUNMatrix J, void *user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    return fJ(t, cj, x, dx, xdot, J, user_data, tmp1, tmp2, tmp3);
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
int IDASolver::fJBandB(realtype t, realtype cj, N_Vector x, N_Vector dx,
                       N_Vector xB, N_Vector dxB, N_Vector xBdot, SUNMatrix JB,
                       void *user_data, N_Vector tmp1B, N_Vector tmp2B,
                       N_Vector tmp3B) {
    return fJB(t, cj, x, dx, xB, dxB, xBdot, JB, user_data, tmp1B, tmp2B,
               tmp3B);
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
    int IDASolver::fJv(realtype t, N_Vector x, N_Vector dx, N_Vector  /*xdot*/, N_Vector v, N_Vector Jv,
                   realtype cj, void *user_data, N_Vector  /*tmp1*/, N_Vector  /*tmp2*/) {
        auto model = static_cast<Model_DAE*>(user_data);
        model->fJv(t, x, dx, v, Jv, cj);
        return model->checkFinite(model->nx_solver,N_VGetArrayPointer(Jv),"Jacobian");
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
    int IDASolver::fJvB(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector  /*xBdot*/,
                    N_Vector vB, N_Vector JvB, realtype cj, void *user_data,
                    N_Vector  /*tmpB1*/, N_Vector  /*tmpB2*/) {
        auto model = static_cast<Model_DAE*>(user_data);
        model->fJvB(t, x, dx, xB, dxB, vB, JvB, cj);
        return model->checkFinite(model->nx_solver,N_VGetArrayPointer(JvB),"Jacobian");
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
        auto model = static_cast<Model_DAE*>(user_data);
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
        auto model = static_cast<Model_DAE *>(user_data);

        if (t > 1e200 && !amici::checkFinite(model->nx_solver,
                                             N_VGetArrayPointer(x), "fxdot"))
            return AMICI_UNRECOVERABLE_ERROR;
        /* when t is large (typically ~1e300), CVODES may pass all NaN x
           to fxdot from which we typically cannot recover. To save time
           on normal execution, we do not always want to check finiteness
           of x, but only do so when t is large and we expect problems. */

        model->fxdot(t, x, dx, xdot);
        return model->checkFinite(model->nx_solver, N_VGetArrayPointer(xdot),
                                  "fxdot");
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
        auto model = static_cast<Model_DAE*>(user_data);
        model->fxBdot(t, x, dx, xB, dxB, xBdot);
        return model->checkFinite(model->nx_solver,N_VGetArrayPointer(xBdot),"xBdot");
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
        auto model = static_cast<Model_DAE*>(user_data);
        model->fqBdot(t, x, dx, xB, dxB, qBdot);
        return model->checkFinite(model->nJ*model->nplist(),N_VGetArrayPointer(qBdot),"qBdot");
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
    int IDASolver::fsxdot(int  /*Ns*/, realtype t, N_Vector x, N_Vector dx, N_Vector  /*xdot*/,
                      N_Vector *sx, N_Vector *sdx, N_Vector *sxdot, void *user_data,
                      N_Vector  /*tmp1*/, N_Vector  /*tmp2*/, N_Vector  /*tmp3*/) {
        auto model = static_cast<Model_DAE*>(user_data);
        for(int ip = 0; ip < model->nplist(); ip++){
            model->fsxdot(t, x, dx, ip, sx[ip], sdx[ip], sxdot[ip]);
            if(model->checkFinite(model->nx_solver,N_VGetArrayPointer(sxdot[ip]),"sxdot") != AMICI_SUCCESS)
                return AMICI_RECOVERABLE_ERROR;
        }
        return AMICI_SUCCESS;
    }

} // namespace amici
