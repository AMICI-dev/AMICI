#include "amici/misc.h"
#include "amici/model_ode.h"
#include "amici/solver_cvodes.h"
#include "amici/exception.h"

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_impl.h>
#include <cvodes/cvodes_diag.h>

#include "amici/sundials_linsol_wrapper.h"


#include <amd.h>
#include <btf.h>
#include <colamd.h>
#include <klu.h>

#define ZERO    RCONST(0.0)
#define ONE     RCONST(1.0)
#define FOUR    RCONST(4.0)

/**
 * @ brief extract information from a property of a matlab class (matrix)
 * @ param OPTION name of the property
 */

namespace amici {

// Ensure AMICI options are in sync with Sundials options
static_assert((int)InternalSensitivityMethod::simultaneous == CV_SIMULTANEOUS, "");
static_assert((int)InternalSensitivityMethod::staggered == CV_STAGGERED, "");
static_assert((int)InternalSensitivityMethod::staggered1 == CV_STAGGERED1, "");

static_assert((int)InterpolationType::hermite == CV_HERMITE, "");
static_assert((int)InterpolationType::polynomial == CV_POLYNOMIAL, "");

static_assert((int)LinearMultistepMethod::adams == CV_ADAMS, "");
static_assert((int)LinearMultistepMethod::BDF == CV_BDF, "");

static_assert(AMICI_ROOT_RETURN == CV_ROOT_RETURN, "");


void CVodeSolver::init(AmiVector *x, AmiVector * /*dx*/, realtype t) {
    int status = CVodeInit(solverMemory.get(), fxdot, t, x->getNVector());
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeInit");
}

void CVodeSolver::binit(int which, AmiVector *xB, AmiVector * /*dxB*/, realtype t) {
    int status = CVodeInitB(solverMemory.get(), which, fxBdot, t, xB->getNVector());
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeInitB");
}

void CVodeSolver::qbinit(int which, AmiVector *qBdot) {
    int status = CVodeQuadInitB(solverMemory.get(), which, fqBdot, qBdot->getNVector());
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeQuadInitB");
}

void CVodeSolver::rootInit(int ne) {
    int status = CVodeRootInit(solverMemory.get(), ne, froot);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeRootInit");
}

void CVodeSolver::sensInit1(AmiVectorArray *sx, AmiVectorArray * /*sdx*/, int nplist) {
    int status = CVodeSensInit1(solverMemory.get(), nplist, static_cast<int>(getSensitivityMethod()), fsxdot,
                          sx->getNVectorArray());
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSensInit1");
}

void CVodeSolver::setDenseJacFn() {
    int status = CVodeSetJacFn(solverMemory.get(), fJ);
    if(status != CV_SUCCESS)
        throw CvodeException(status,"CVodeSetJacFn");
}

void CVodeSolver::setSparseJacFn() {
    int status = CVodeSetJacFn(solverMemory.get(), fJSparse);
    if(status != CV_SUCCESS)
        throw CvodeException(status,"CVodeSetJacFn");
}

void CVodeSolver::setBandJacFn() {
    int status = CVodeSetJacFn(solverMemory.get(), fJBand);
    if(status != CV_SUCCESS)
        throw CvodeException(status,"CVodeSetJacFn");
}

void CVodeSolver::setJacTimesVecFn() {
    int status = CVodeSetJacTimes(solverMemory.get(), nullptr, fJv);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetJacTimes");
}

void CVodeSolver::setDenseJacFnB(int which) {
    int status = CVodeSetJacFnB(solverMemory.get(), which, fJB);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetJacFnB");
}

void CVodeSolver::setSparseJacFnB(int which) {
    int status = CVodeSetJacFnB(solverMemory.get(), which, fJSparseB);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetJacFnB");
}

void CVodeSolver::setBandJacFnB(int which) {
    int status = CVodeSetJacFnB(solverMemory.get(), which, fJBandB);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetJacFnB");
}

void CVodeSolver::setJacTimesVecFnB(int which) {
    int status = CVodeSetJacTimesB(solverMemory.get(), which, nullptr, fJvB);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetJacTimesB");
}

Solver *CVodeSolver::clone() const {
    return new CVodeSolver(*this);
}

void CVodeSolver::allocateSolver() {
    solverMemory = std::unique_ptr<void, std::function<void(void *)>>
    (CVodeCreate(static_cast<int>(lmm)),
                   [](void *ptr) { CVodeFree(&ptr); });
}

void CVodeSolver::setSStolerances(double rtol, double atol) {
    int status = CVodeSStolerances(solverMemory.get(), rtol, atol);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSStolerances");
}

void CVodeSolver::setSensSStolerances(double rtol, double *atol) {
    int status = CVodeSensSStolerances(solverMemory.get(), rtol, atol);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSensEEtolerances");
}

void CVodeSolver::setSensErrCon(bool error_corr) {
    int status = CVodeSetSensErrCon(solverMemory.get(), error_corr);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetSensErrCon");
}

void CVodeSolver::setQuadErrConB(int which, bool flag) {
    int status = CVodeSetQuadErrConB(solverMemory.get(), which, flag);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetQuadErrConB");
}

void CVodeSolver::getRootInfo(int *rootsfound) const {
    int status = CVodeGetRootInfo(solverMemory.get(), rootsfound);
    if(status != CV_SUCCESS)
        throw CvodeException(status,"CVodeGetRootInfo");
}

void CVodeSolver::setLinearSolver()
{
    int status = CVodeSetLinearSolver(solverMemory.get(), linearSolver->get(), linearSolver->getMatrix());
    if(status != CV_SUCCESS)
        throw CvodeException(status,"setLinearSolver");
}

void CVodeSolver::setLinearSolverB(int which)
{
    int status = CVodeSetLinearSolverB(solverMemory.get(), which, linearSolverB->get(), linearSolverB->getMatrix());
    if(status != CV_SUCCESS)
        throw CvodeException(status,"setLinearSolverB");
}

void CVodeSolver::setNonLinearSolver()
{
    int status = CVodeSetNonlinearSolver(solverMemory.get(), nonLinearSolver->get());
    if(status != CV_SUCCESS)
        throw CvodeException(status,"CVodeSetNonlinearSolver");
}

void CVodeSolver::setNonLinearSolverSens()
{
    if(getSensitivityOrder() < SensitivityOrder::first)
        return;
    if(getSensitivityMethod() != SensitivityMethod::forward)
        return;

    int status = CV_SUCCESS;

    switch (ism) {
    case InternalSensitivityMethod::staggered:
        status = CVodeSetNonlinearSolverSensStg(solverMemory.get(), nonLinearSolverSens->get());
        break;
    case InternalSensitivityMethod::simultaneous:
        status = CVodeSetNonlinearSolverSensSim(solverMemory.get(), nonLinearSolverSens->get());
        break;
    case InternalSensitivityMethod::staggered1:
        status = CVodeSetNonlinearSolverSensStg1(solverMemory.get(), nonLinearSolverSens->get());
        break;
    default:
        throw AmiException("Unsupported internal sensitivity method selected: %d", ism);
    }

    if(status != CV_SUCCESS)
        throw CvodeException(status,"CVodeSolver::setNonLinearSolverSens");

}

void CVodeSolver::setNonLinearSolverB(int which)
{
    int status = CVodeSetNonlinearSolverB(solverMemory.get(), which, nonLinearSolverB->get());
    if(status != CV_SUCCESS)
        throw CvodeException(status,"CVodeSetNonlinearSolverB");
}

void CVodeSolver::setErrHandlerFn() {
    int status = CVodeSetErrHandlerFn(solverMemory.get(), wrapErrHandlerFn, nullptr);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetErrHandlerFn");
}

void CVodeSolver::setUserData(Model *model) {
    int status = CVodeSetUserData(solverMemory.get(), model);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetUserData");
}

void CVodeSolver::setUserDataB(int which, Model *model) {
    int status = CVodeSetUserDataB(solverMemory.get(), which, model);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetUserDataB");
}

void CVodeSolver::setMaxNumSteps(long mxsteps) {
    int status = CVodeSetMaxNumSteps(solverMemory.get(), mxsteps);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetMaxNumSteps");
}

void CVodeSolver::setStabLimDet(int stldet) {
    int status = CVodeSetStabLimDet(solverMemory.get(), stldet);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetStabLimDet");
}

void CVodeSolver::setStabLimDetB(int which, int stldet) {
    int status = CVodeSetStabLimDetB(solverMemory.get(), which, stldet);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetStabLimDetB");
}

void CVodeSolver::setId(Model *model) { }

void CVodeSolver::setSuppressAlg(bool flag) { }

void CVodeSolver::resetState(void *ami_mem, N_Vector y0) {

    auto cv_mem = static_cast<CVodeMem>(ami_mem);
    /* here we force the order in the next step to zero, and update the
     Nordsieck history array, this is largely copied from CVodeReInit with
     explanations from cvodes_impl.h
     */

    /* Set step parameters */

    /* current order */
    cv_mem->cv_q      = 1;
    /* L = q + 1 */
    cv_mem->cv_L      = 2;
    /* number of steps to wait before updating in q */
    cv_mem->cv_qwait  = cv_mem->cv_L;
    /* last successful q value used                */
    cv_mem->cv_qu     = 0;
    /* last successful h value used                */
    cv_mem->cv_hu     = ZERO;
    /* tolerance scale factor                      */
    cv_mem->cv_tolsf  = ONE;

    /* Initialize other integrator optional outputs */

    /* actual initial stepsize                     */
    cv_mem->cv_h0u      = ZERO;
    /* step size to be used on the next step       */
    cv_mem->cv_next_h   = ZERO;
    /* order to be used on the next step           */
    cv_mem->cv_next_q   = 0;

    /* write updated state to Nordsieck history array  */
    N_VScale(ONE, y0, cv_mem->cv_zn[0]);
}


void CVodeSolver::reInitPostProcessF(realtype *t, AmiVector *yout,
                                     AmiVector * /*ypout*/, realtype tnext) {
    reInitPostProcess(solverMemory.get(), t, yout, tnext);
}

void CVodeSolver::reInitPostProcessB(int which, realtype *t, AmiVector *yBout,
                                     AmiVector * /*ypBout*/, realtype tnext) {
    reInitPostProcess(CVodeGetAdjCVodeBmem(solverMemory.get(), which),
                      t, yBout, tnext);
}

void CVodeSolver::reInitPostProcess(void *ami_mem, realtype *t,
                                    AmiVector *yout, realtype tout) {
    auto cv_mem = static_cast<CVodeMem>(ami_mem);
    auto nst_tmp = cv_mem->cv_nst;
    cv_mem->cv_nst = 0;

    auto status = CVodeSetStopTime(cv_mem, tout);
    if(status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetStopTime");

    status = CVode(ami_mem, tout, yout->getNVector(),
                        t, CV_ONE_STEP);

    if(status != CV_SUCCESS)
        throw CvodeException(status, "reInitPostProcess");

    cv_mem->cv_nst = nst_tmp+1;
    if (cv_mem->cv_adjMallocDone == SUNTRUE) {
        /* add new step to history array, this is copied from CVodeF */
        auto ca_mem = cv_mem->cv_adj_mem;
        auto dt_mem = ca_mem->dt_mem;

        if (cv_mem->cv_nst % ca_mem->ca_nsteps == 0) {
            /* currently not implemented, we should never get here as we
             limit cv_mem->cv_nst < ca_mem->ca_nsteps, keeping this for
             future regression */
            throw CvodeException(AMICI_ERROR, "reInitPostProcess");
        }

        /* Load next point in dt_mem */
        dt_mem[cv_mem->cv_nst % ca_mem->ca_nsteps]->t = *t;
        ca_mem->ca_IMstore(cv_mem, dt_mem[cv_mem->cv_nst % ca_mem->ca_nsteps]);

        /* Set t1 field of the current ckeck point structure
         for the case in which there will be no future
         check points */
        ca_mem->ck_mem->ck_t1 = *t;

        /* tfinal is now set to *tret */
        ca_mem->ca_tfinal = *t;
    }
}

void CVodeSolver::reInit(realtype t0, AmiVector *yy0, AmiVector * /*yp0*/) {
    auto cv_mem = static_cast<CVodeMem>(solverMemory.get());
    /* set time */
    cv_mem->cv_tn = t0;
    resetState(cv_mem, yy0->getNVector());
}

void CVodeSolver::sensReInit(AmiVectorArray *yS0, AmiVectorArray * /*ypS0*/) {
    auto cv_mem = static_cast<CVodeMem>(solverMemory.get());
    /* Initialize znS[0] in the history array */

    for (int is=0; is<cv_mem->cv_Ns; is++)
        cv_mem->cv_cvals[is] = ONE;

    int status = N_VScaleVectorArray(cv_mem->cv_Ns, cv_mem->cv_cvals,
                                 yS0->getNVectorArray(), cv_mem->cv_znS[0]);
    if(status != CV_SUCCESS)
         throw CvodeException(CV_VECTOROP_ERR,"CVodeSensReInit");
}

void CVodeSolver::setSensParams(realtype *p, realtype *pbar, int *plist) {
    int status = CVodeSetSensParams(solverMemory.get(), p, pbar, plist);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetSensParams");
}

void CVodeSolver::getDky(realtype t, int k, AmiVector *dky) const {
    int status = CVodeGetDky(solverMemory.get(), t, k, dky->getNVector());
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeGetDky");
}

void CVodeSolver::getSens(realtype *tret, AmiVectorArray *yySout) const {
    int status = CVodeGetSens(solverMemory.get(), tret, yySout->getNVectorArray());
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeGetSens");
}

void CVodeSolver::adjInit() {
    int status = CVodeAdjInit(solverMemory.get(), static_cast<int>(maxsteps), static_cast<int>(interpType));
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeAdjInit");
}

void CVodeSolver::allocateSolverB(int *which) {
    int status = CVodeCreateB(solverMemory.get(), static_cast<int>(lmm), which);

    if (*which + 1 > static_cast<int>(solverMemoryB.size()))
        solverMemoryB.resize(*which + 1);
    solverMemoryB.at(*which) = std::unique_ptr<void, std::function<void(void *)>>
    (getAdjBmem(solverMemory.get(), *which), [](void *ptr){});
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeCreateB");
}

void CVodeSolver::reInitB(int which, realtype tB0, AmiVector *yyB0,
                            AmiVector * /*ypB0*/) {
    auto cv_memB = static_cast<CVodeMem>(CVodeGetAdjCVodeBmem(solverMemory.get(), which));
    cv_memB->cv_tn = tB0;
    resetState(cv_memB, yyB0->getNVector());
}

void CVodeSolver::setSStolerancesB(int which, realtype relTolB,
                                  realtype absTolB) {
    int status = CVodeSStolerancesB(solverMemory.get(), which, relTolB, absTolB);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSStolerancesB");
}

void CVodeSolver::quadReInitB(int which, AmiVector *yQB0) {
    auto cv_memB = static_cast<CVodeMem>(CVodeGetAdjCVodeBmem(solverMemory.get(), which));
    N_VScale(ONE, yQB0->getNVector(), cv_memB->cv_znQ[0]);
}

void CVodeSolver::quadSStolerancesB(int which, realtype reltolQB,
                                      realtype abstolQB) {
    int status = CVodeQuadSStolerancesB(solverMemory.get(), which, reltolQB, abstolQB);
    if(status != CV_SUCCESS)
        throw CvodeException(status,"CVodeQuadSStolerancesB");
}

void CVodeSolver::getQuadB(int which, realtype *tret, AmiVector *qB) const {
   int status = CVodeGetQuadB(solverMemory.get(), which, tret, qB->getNVector());
   if(status != CV_SUCCESS)
        throw CvodeException(status,"CVodeGetQuadB");
}


int CVodeSolver::solve(realtype tout, AmiVector *yret, AmiVector * /*ypret*/,
                          realtype *tret, int itask) {
    int status = CVode(solverMemory.get(), tout, yret->getNVector(), tret, itask);
    if(status<0) {
        throw IntegrationFailure(status,*tret);
    }

    solverWasCalled = true;
    return status;
}

int CVodeSolver::solveF(realtype tout, AmiVector *yret, AmiVector * /*ypret*/,
                           realtype *tret, int itask, int *ncheckPtr) {
    int status = CVodeF(solverMemory.get(), tout, yret->getNVector(), tret, itask, ncheckPtr);
    if(status<0) {
        throw IntegrationFailure(status,*tret);
    }

    solverWasCalled = true;
    return status;
}

void CVodeSolver::solveB(realtype tBout, int itaskB) {
    int status = CVodeB(solverMemory.get(), tBout, itaskB);
    if(status != CV_SUCCESS)
         throw IntegrationFailureB(status,tBout);
}

void CVodeSolver::setMaxNumStepsB(int which, long mxstepsB) {
    int status = CVodeSetMaxNumStepsB(solverMemory.get(), which, mxstepsB);
    if(status != CV_SUCCESS)
        throw CvodeException(status,"CVodeSetMaxNumStepsB");
}

void CVodeSolver::diag()
{
    int status = CVDiag(solverMemory.get());
    if(status != CV_SUCCESS)
        throw CvodeException(status,"CVDiag");
}

void CVodeSolver::diagB(int which)
{
    int status = CVDiagB(solverMemory.get(), which);
    if(status != CV_SUCCESS)
        throw CvodeException(status,"CVDiagB");

}

void CVodeSolver::getB(int which, realtype *tret, AmiVector *yy, AmiVector * /*yp*/) const {
    int status = CVodeGetB(solverMemory.get(), which, tret, yy->getNVector());
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeGetB");
}


void CVodeSolver::getNumSteps(void *ami_mem, long *numsteps) const {
    int status = CVodeGetNumSteps(ami_mem, numsteps);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeGetNumSteps");
}

void CVodeSolver::getNumRhsEvals(void *ami_mem, long *numrhsevals) const {
    int status = CVodeGetNumRhsEvals(ami_mem, numrhsevals);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeGetNumRhsEvals");
}

void CVodeSolver::getNumErrTestFails(void *ami_mem, long *numerrtestfails) const {
    int status = CVodeGetNumErrTestFails(ami_mem, numerrtestfails);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeGetNumErrTestFails");
}

void CVodeSolver::getNumNonlinSolvConvFails(void *ami_mem,
                                              long *numnonlinsolvconvfails) const {
    int status = CVodeGetNumNonlinSolvConvFails(ami_mem, numnonlinsolvconvfails);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeGetNumNonlinSolvConvFails");
}

void CVodeSolver::getLastOrder(void *ami_mem, int *order) const {
    int status = CVodeGetLastOrder(ami_mem, order);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeGetLastOrder");
}

void *CVodeSolver::getAdjBmem(void *ami_mem, int which) {
    return CVodeGetAdjCVodeBmem(ami_mem, which);
}

void CVodeSolver::calcIC(realtype tout1, AmiVector *x, AmiVector *dx) { };

void CVodeSolver::calcICB(int which, realtype tout1, AmiVector *xB,
                             AmiVector *dxB) {};

void CVodeSolver::setStopTime(realtype tstop) {
    int status = CVodeSetStopTime(solverMemory.get(), tstop);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeSetStopTime");
}

void CVodeSolver::turnOffRootFinding() {
    int status = CVodeRootInit(solverMemory.get(), 0, nullptr);
    if(status != CV_SUCCESS)
         throw CvodeException(status,"CVodeRootInit");
}

int CVodeSolver::nplist() const {
    if (!solverMemory)
        throw AmiException("Solver has not been allocated, information is not available");
    auto cv_mem = static_cast<CVodeMem>(solverMemory.get());
    return cv_mem->cv_Ns;
}

int CVodeSolver::nx() const {
    if (!solverMemory)
        throw AmiException("Solver has not been allocated, information is not available");
    auto cv_mem = static_cast<CVodeMem>(solverMemory.get());
    return NV_LENGTH_S(cv_mem->cv_zn[0]);
}

const Model *CVodeSolver::getModel() const {
    if (!solverMemory)
        throw AmiException("Solver has not been allocated, information is not available");
    auto cv_mem = static_cast<CVodeMem>(solverMemory.get());
    return static_cast<Model *>(cv_mem->cv_user_data);
}

bool CVodeSolver::getMallocDone() const {
    if (!solverMemory)
        return false;
    auto cv_mem = static_cast<CVodeMem>(solverMemory.get());
    return cv_mem->cv_MallocDone;
}

bool CVodeSolver::getAdjMallocDone() const {
    if (!solverMemory)
        return false;
    auto cv_mem = static_cast<CVodeMem>(solverMemory.get());
    return cv_mem->cv_adjMallocDone;
}

/** Jacobian of xdot with respect to states x
 * @param N number of state variables
 * @param t timepoint
 * @param x Vector with the states
 * @param xdot Vector with the right hand side
 * @param J Matrix to which the Jacobian will be written
 * @param user_data object with user input @type Model_ODE
 * @param tmp1 temporary storage vector
 * @param tmp2 temporary storage vector
 * @param tmp3 temporary storage vector
 * @return status flag indicating successful execution
 **/
int CVodeSolver::fJ(realtype t, N_Vector x, N_Vector xdot, SUNMatrix J,
                    void *user_data, N_Vector  /*tmp1*/, N_Vector  /*tmp2*/,
                    N_Vector  /*tmp3*/) {
    auto model = static_cast<Model_ODE *>(user_data);
    model->fJ(t, x, xdot, J);
    return model->checkFinite(SM_ROWS_D(J), SM_DATA_D(J), "Jacobian");
}

/** Jacobian of xBdot with respect to adjoint state xB
 * @param NeqBdot number of adjoint state variables
 * @param t timepoint
 * @param x Vector with the states
 * @param xB Vector with the adjoint states
 * @param xBdot Vector with the adjoint right hand side
 * @param JB Matrix to which the Jacobian will be written
 * @param user_data object with user input @type Model_ODE
 * @param tmp1B temporary storage vector
 * @param tmp2B temporary storage vector
 * @param tmp3B temporary storage vector
 * @return status flag indicating successful execution
 **/
int CVodeSolver::fJB(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot,
                     SUNMatrix JB, void *user_data, N_Vector  /*tmp1B*/,
                     N_Vector  /*tmp2B*/, N_Vector  /*tmp3B*/) {
    auto model = static_cast<Model_ODE *>(user_data);
    model->fJB(t, x, xB, xBdot, JB);
    return model->checkFinite(SM_ROWS_D(JB), SM_DATA_D(JB), "Jacobian");
}

/** J in sparse form (for sparse solvers from the SuiteSparse Package)
 * @param t timepoint
 * @param x Vector with the states
 * @param xdot Vector with the right hand side
 * @param J Matrix to which the Jacobian will be written
 * @param user_data object with user input @type Model_ODE
 * @param tmp1 temporary storage vector
 * @param tmp2 temporary storage vector
 * @param tmp3 temporary storage vector
 * @return status flag indicating successful execution
 */
int CVodeSolver::fJSparse(realtype t, N_Vector x, N_Vector  /*xdot*/, SUNMatrix J,
                          void *user_data, N_Vector  /*tmp1*/, N_Vector  /*tmp2*/,
                          N_Vector  /*tmp3*/) {
    auto model = static_cast<Model_ODE *>(user_data);
    model->fJSparse(t, x, J);
    return model->checkFinite(SM_NNZ_S(J), SM_DATA_S(J), "Jacobian");
}

/** JB in sparse form (for sparse solvers from the SuiteSparse Package)
 * @param t timepoint
 * @param x Vector with the states
 * @param xB Vector with the adjoint states
 * @param xBdot Vector with the adjoint right hand side
 * @param JB Matrix to which the Jacobian will be written
 * @param user_data object with user input @type Model_ODE
 * @param tmp1B temporary storage vector
 * @param tmp2B temporary storage vector
 * @param tmp3B temporary storage vector
 * @return status flag indicating successful execution
 */
int CVodeSolver::fJSparseB(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot,
                           SUNMatrix JB, void *user_data, N_Vector  /*tmp1B*/,
                           N_Vector  /*tmp2B*/, N_Vector  /*tmp3B*/) {
    auto model = static_cast<Model_ODE *>(user_data);
    model->fJSparseB(t, x, xB, xBdot, JB);
    return model->checkFinite(SM_NNZ_S(JB), SM_DATA_S(JB), "Jacobian");
}

/** J in banded form (for banded solvers)
 * @param N number of states
 * @param mupper upper matrix bandwidth
 * @param mlower lower matrix bandwidth
 * @param t timepoint
 * @param x Vector with the states
 * @param xdot Vector with the right hand side
 * @param J Matrix to which the Jacobian will be written
 * @param user_data object with user input @type Model_ODE
 * @param tmp1 temporary storage vector
 * @param tmp2 temporary storage vector
 * @param tmp3 temporary storage vector
 * @return status flag indicating successful execution
 */
int CVodeSolver::fJBand(realtype t, N_Vector x, N_Vector xdot, SUNMatrix J,
                        void *user_data, N_Vector tmp1, N_Vector tmp2,
                        N_Vector tmp3) {
    return fJ(t, x, xdot, J, user_data, tmp1, tmp2, tmp3);
}

/** JB in banded form (for banded solvers)
 * @param NeqBdot number of states
 * @param mupper upper matrix bandwidth
 * @param mlower lower matrix bandwidth
 * @param t timepoint
 * @param x Vector with the states
 * @param xB Vector with the adjoint states
 * @param xBdot Vector with the adjoint right hand side
 * @param JB Matrix to which the Jacobian will be written
 * @param user_data object with user input @type Model_ODE
 * @param tmp1B temporary storage vector
 * @param tmp2B temporary storage vector
 * @param tmp3B temporary storage vector
 * @return status flag indicating successful execution
 */
int CVodeSolver::fJBandB(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot,
                         SUNMatrix JB, void *user_data, N_Vector tmp1B,
                         N_Vector tmp2B, N_Vector tmp3B) {
    return fJB(t, x, xB, xBdot, JB, user_data, tmp1B, tmp2B, tmp3B);
}

/** diagonalized Jacobian (for preconditioning)
 * @param t timepoint
 * @param JDiag Vector to which the Jacobian diagonal will be written
 * @param x Vector with the states
 * @param user_data object with user input @type Model_ODE
 **/
int CVodeSolver::fJDiag(realtype t, N_Vector JDiag, N_Vector x,
                        void *user_data) {
    auto model = static_cast<Model_ODE *>(user_data);
    model->fJDiag(t, JDiag, x);
    return model->checkFinite(model->nx_solver, N_VGetArrayPointer(JDiag),
                              "Jacobian");
}

    /** Matrix vector product of J with a vector v (for iterative solvers)
     * @param t timepoint
     * @param x Vector with the states
     * @param xdot Vector with the right hand side
     * @param v Vector with which the Jacobian is multiplied
     * @param Jv Vector to which the Jacobian vector product will be
     *written
     * @param user_data object with user input @type Model_ODE
     * @param tmp temporary storage vector
     * @return status flag indicating successful execution
     **/
    int CVodeSolver::fJv(N_Vector v, N_Vector Jv, realtype t, N_Vector x, N_Vector  /*xdot*/,
                   void *user_data, N_Vector  /*tmp*/) {
        auto model = static_cast<Model_ODE*>(user_data);
        model->fJv(v,Jv,t,x);
        return model->checkFinite(model->nx_solver,N_VGetArrayPointer(Jv),"Jacobian");
    }

    /** Matrix vector product of JB with a vector v (for iterative solvers)
     * @param t timepoint
     * @param x Vector with the states
     * @param xB Vector with the adjoint states
     * @param xBdot Vector with the adjoint right hand side
     * @param vB Vector with which the Jacobian is multiplied
     * @param JvB Vector to which the Jacobian vector product will be
     *written
     * @param user_data object with user input @type Model_ODE
     * @param tmpB temporary storage vector
     * @return status flag indicating successful execution
     **/
    int CVodeSolver::fJvB(N_Vector vB, N_Vector JvB, realtype t, N_Vector x, N_Vector xB, N_Vector  /*xBdot*/,
                    void *user_data, N_Vector  /*tmpB*/) {
        auto model = static_cast<Model_ODE*>(user_data);
        model->fJvB(vB, JvB, t, x, xB);
        return model->checkFinite(model->nx_solver,N_VGetArrayPointer(JvB),"Jacobian");
    }

    /** Event trigger function for events
     * @param t timepoint
     * @param x Vector with the states
     * @param root array with root function values
     * @param user_data object with user input @type Model_ODE
     * @return status flag indicating successful execution
     */
    int CVodeSolver::froot(realtype t, N_Vector x, realtype *root,
                     void *user_data) {
        auto model = static_cast<Model_ODE*>(user_data);
        model->froot(t, x, root);
        return model->checkFinite(model->ne,root,"root function");
    }

    /** residual function of the ODE
     * @param t timepoint
     * @param x Vector with the states
     * @param xdot Vector with the right hand side
     * @param user_data object with user input @type Model_ODE
     * @return status flag indicating successful execution
     */
    int CVodeSolver::fxdot(realtype t, N_Vector x, N_Vector xdot,
                           void *user_data) {
        auto model = static_cast<Model_ODE *>(user_data);

        if (t > 1e200 && !amici::checkFinite(model->nx_solver,
                                             N_VGetArrayPointer(x), "fxdot"))
            return AMICI_UNRECOVERABLE_ERROR;
            /* when t is large (typically ~1e300), CVODES may pass all NaN x
               to fxdot from which we typically cannot recover. To save time
               on normal execution, we do not always want to check finiteness
               of x, but only do so when t is large and we expect problems. */

        model->fxdot(t, x, xdot);
        return model->checkFinite(model->nx_solver, N_VGetArrayPointer(xdot),
                                  "fxdot");
    }

    /** Right hand side of differential equation for adjoint state xB
     * @param t timepoint
     * @param x Vector with the states
     * @param xB Vector with the adjoint states
     * @param xBdot Vector with the adjoint right hand side
     * @param user_data object with user input @type Model_ODE
     * @return status flag indicating successful execution
     */
    int CVodeSolver::fxBdot(realtype t, N_Vector x, N_Vector xB,
                      N_Vector xBdot, void *user_data) {
        auto model = static_cast<Model_ODE*>(user_data);
        model->fxBdot(t, x, xB, xBdot);
        return model->checkFinite(model->nx_solver,N_VGetArrayPointer(xBdot),"fxBdot");
    }

    /** Right hand side of integral equation for quadrature states qB
     * @param t timepoint
     * @param x Vector with the states
     * @param xB Vector with the adjoint states
     * @param qBdot Vector with the adjoint quadrature right hand side
     * @param user_data pointer to temp data object
     * @return status flag indicating successful execution
     */
    int CVodeSolver::fqBdot(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot,
                      void *user_data) {
        auto model = static_cast<Model_ODE*>(user_data);
        model->fqBdot(t, x, xB, qBdot);
        return model->checkFinite(model->nplist()*model->nJ,N_VGetArrayPointer(qBdot),"qBdot");
    }

    /** Right hand side of differential equation for state sensitivities sx
     * @param Ns number of parameters
     * @param t timepoint
     * @param x Vector with the states
     * @param xdot Vector with the right hand side
     * @param ip parameter index
     * @param sx Vector with the state sensitivities
     * @param sxdot Vector with the sensitivity right hand side
     * @param user_data object with user input @type Model_ODE
     * @param tmp1 temporary storage vector
     * @param tmp2 temporary storage vector
     * @param tmp3 temporary storage vector
     * @return status flag indicating successful execution
     */
    int CVodeSolver::fsxdot(int  /*Ns*/, realtype t, N_Vector x, N_Vector  /*xdot*/, int ip,
                      N_Vector sx, N_Vector sxdot, void *user_data,
                      N_Vector  /*tmp1*/, N_Vector  /*tmp2*/) {
        auto model = static_cast<Model_ODE*>(user_data);
        model->fsxdot(t, x, ip, sx, sxdot);
        return model->checkFinite(model->nx_solver,N_VGetArrayPointer(sxdot),"sxdot");
    }

    bool operator ==(const CVodeSolver &a, const CVodeSolver &b)
    {
        return static_cast<Solver const&>(a) == static_cast<Solver const&>(b);
    }
} // namespace amici
