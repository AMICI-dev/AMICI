#include "amici/solver_cvodes.h"

#include "amici/exception.h"
#include "amici/misc.h"
#include "amici/model_ode.h"
#include "amici/sundials_linsol_wrapper.h"

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_diag.h>
#include <cvodes/cvodes_impl.h>

#include <amd.h>
#include <btf.h>
#include <colamd.h>
#include <klu.h>

#define ZERO RCONST(0.0)
#define ONE RCONST(1.0)
#define FOUR RCONST(4.0)

namespace amici {

// Ensure AMICI options are in sync with Sundials options
static_assert((int)InternalSensitivityMethod::simultaneous == CV_SIMULTANEOUS,
              "");
static_assert((int)InternalSensitivityMethod::staggered == CV_STAGGERED, "");
static_assert((int)InternalSensitivityMethod::staggered1 == CV_STAGGERED1, "");

static_assert((int)InterpolationType::hermite == CV_HERMITE, "");
static_assert((int)InterpolationType::polynomial == CV_POLYNOMIAL, "");

static_assert((int)LinearMultistepMethod::adams == CV_ADAMS, "");
static_assert((int)LinearMultistepMethod::BDF == CV_BDF, "");

static_assert(AMICI_ROOT_RETURN == CV_ROOT_RETURN, "");


/*
 * The following static members are callback function to CVODES.
 * Their signatures must not be changes.
 */
static int fxdot(realtype t, N_Vector x, N_Vector xdot, void *user_data);

static int fJSparse(realtype t, N_Vector x, N_Vector xdot, SUNMatrix J,
                    void *user_data, N_Vector tmp1, N_Vector tmp2,
                    N_Vector tmp3);

static int fJ(realtype t, N_Vector x, N_Vector xdot, SUNMatrix J,
              void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int fJB(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot,
               SUNMatrix JB, void *user_data, N_Vector tmp1B,
               N_Vector tmp2B, N_Vector tmp3B);

static int fJSparseB(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot,
                     SUNMatrix JB, void *user_data, N_Vector tmp1B,
                     N_Vector tmp2B, N_Vector tmp3B);

static int fJBand(realtype t, N_Vector x, N_Vector xdot, SUNMatrix J,
                  void *user_data, N_Vector tmp1, N_Vector tmp2,
                  N_Vector tmp3);

static int fJBandB(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot,
                   SUNMatrix JB, void *user_data, N_Vector tmp1B,
                   N_Vector tmp2B, N_Vector tmp3B);

static int fJv(N_Vector v, N_Vector Jv, realtype t, N_Vector x,
               N_Vector xdot, void *user_data, N_Vector tmp);

static int fJvB(N_Vector vB, N_Vector JvB, realtype t, N_Vector x,
                N_Vector xB, N_Vector xBdot, void *user_data,
                N_Vector tmpB);

static int froot(realtype t, N_Vector x, realtype *root, void *user_data);

static int fxBdot(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot,
                  void *user_data);

static int fqBdot(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot,
                  void *user_data);

static int fsxdot(int Ns, realtype t, N_Vector x, N_Vector xdot, int ip,
                  N_Vector sx, N_Vector sxdot, void *user_data,
                  N_Vector tmp1, N_Vector tmp2);


/* Function implementations */

void CVodeSolver::init(const realtype t0, const AmiVector &x0,
                       const AmiVector & /*dx0*/) const {
    solverWasCalledF = false;
    forceReInitPostProcessF = false;
    t = t0;
    x = x0;
    int status;
    if (getInitDone()) {
        status = CVodeReInit(solverMemory.get(), t0, x.getNVector());
    } else {
        status = CVodeInit(solverMemory.get(), fxdot, t0, x.getNVector());
        setInitDone();
    }
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeInit");
}

void CVodeSolver::sensInit1(const AmiVectorArray &sx0,
                            const AmiVectorArray & /*sdx0*/) const {
    int status = CV_SUCCESS;
    sx = sx0;
    if (getSensitivityMethod() == SensitivityMethod::forward && nplist() > 0) {
        if (getSensInitDone()) {
            status = CVodeSensReInit(
                solverMemory.get(),
                static_cast<int>(getInternalSensitivityMethod()),
                sx.getNVectorArray());
        } else {
            status =
                CVodeSensInit1(solverMemory.get(), nplist(),
                               static_cast<int>(getInternalSensitivityMethod()),
                               fsxdot, sx.getNVectorArray());
            setSensInitDone();
        }
    }
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSensInit1");
}

void CVodeSolver::binit(const int which, const realtype tf,
                        const AmiVector &xB0,
                        const AmiVector & /*dxB0*/) const {
    solverWasCalledB = false;
    forceReInitPostProcessB = false;
    xB = xB0;
    int status;
    if (getInitDoneB(which)) {
        status = CVodeReInitB(solverMemory.get(), which, tf, xB.getNVector());
    } else {
        status =
            CVodeInitB(solverMemory.get(), which, fxBdot, tf, xB.getNVector());
        setInitDoneB(which);
    }
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeInitB");
}

void CVodeSolver::qbinit(const int which, const AmiVector &xQB0) const {
    xQB = xQB0;
    int status;
    if (getQuadInitDoneB(which)) {
        status = CVodeQuadReInitB(solverMemory.get(), which, xQB.getNVector());
    } else {
        status =
            CVodeQuadInitB(solverMemory.get(), which, fqBdot, xQB.getNVector());
        setQuadInitDoneB(which);
    }
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeQuadInitB");
}

void CVodeSolver::rootInit(int ne) const {
    int status = CVodeRootInit(solverMemory.get(), ne, froot);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeRootInit");
}

void CVodeSolver::setDenseJacFn() const {
    int status = CVodeSetJacFn(solverMemory.get(), fJ);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetJacFn");
}

void CVodeSolver::setSparseJacFn() const {
    int status = CVodeSetJacFn(solverMemory.get(), fJSparse);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetJacFn");
}

void CVodeSolver::setBandJacFn() const {
    int status = CVodeSetJacFn(solverMemory.get(), fJBand);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetJacFn");
}

void CVodeSolver::setJacTimesVecFn() const {
    int status = CVodeSetJacTimes(solverMemory.get(), nullptr, fJv);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetJacTimes");
}

void CVodeSolver::setDenseJacFnB(int which) const {
    int status = CVodeSetJacFnB(solverMemory.get(), which, fJB);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetJacFnB");
}

void CVodeSolver::setSparseJacFnB(int which) const {
    int status = CVodeSetJacFnB(solverMemory.get(), which, fJSparseB);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetJacFnB");
}

void CVodeSolver::setBandJacFnB(int which) const {
    int status = CVodeSetJacFnB(solverMemory.get(), which, fJBandB);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetJacFnB");
}

void CVodeSolver::setJacTimesVecFnB(int which) const {
    int status = CVodeSetJacTimesB(solverMemory.get(), which, nullptr, fJvB);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetJacTimesB");
}

Solver *CVodeSolver::clone() const { return new CVodeSolver(*this); }

void CVodeSolver::allocateSolver() const {
    if (!solverMemory)
        solverMemory = std::unique_ptr<void, std::function<void(void *)>>(
            CVodeCreate(static_cast<int>(lmm)),
            [](void *ptr) { CVodeFree(&ptr); });
}

void CVodeSolver::setSStolerances(const double rtol, const double atol) const {
    int status = CVodeSStolerances(solverMemory.get(), rtol, atol);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSStolerances");
}

void CVodeSolver::setSensSStolerances(const double rtol,
                                      const double *atol) const {
    int status = CVodeSensSStolerances(solverMemory.get(), rtol,
                                       const_cast<double *>(atol));
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSensEEtolerances");
}

void CVodeSolver::setSensErrCon(const bool error_corr) const {
    int status = CVodeSetSensErrCon(solverMemory.get(), error_corr);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetSensErrCon");
}

void CVodeSolver::setQuadErrConB(const int which, const bool flag) const {
    int status = CVodeSetQuadErrConB(solverMemory.get(), which, flag);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetQuadErrConB");
}

void CVodeSolver::setQuadErrCon(const bool flag) const {
    int status = CVodeSetQuadErrCon(solverMemory.get(), flag);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetQuadErrCon");
}

void CVodeSolver::getRootInfo(int *rootsfound) const {
    int status = CVodeGetRootInfo(solverMemory.get(), rootsfound);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetRootInfo");
}

void CVodeSolver::setLinearSolver() const {
    int status = CVodeSetLinearSolver(solverMemory.get(), linearSolver->get(),
                                      linearSolver->getMatrix());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "setLinearSolver");
}

void CVodeSolver::setLinearSolverB(int which) const {
    int status =
        CVodeSetLinearSolverB(solverMemory.get(), which, linearSolverB->get(),
                              linearSolverB->getMatrix());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "setLinearSolverB");
}

void CVodeSolver::setNonLinearSolver() const {
    int status =
        CVodeSetNonlinearSolver(solverMemory.get(), nonLinearSolver->get());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetNonlinearSolver");
}

void CVodeSolver::setNonLinearSolverSens() const {
    if (getSensitivityOrder() < SensitivityOrder::first)
        return;
    if (getSensitivityMethod() != SensitivityMethod::forward)
        return;

    int status = CV_SUCCESS;

    switch (ism) {
    case InternalSensitivityMethod::staggered:
        status = CVodeSetNonlinearSolverSensStg(solverMemory.get(),
                                                nonLinearSolverSens->get());
        break;
    case InternalSensitivityMethod::simultaneous:
        status = CVodeSetNonlinearSolverSensSim(solverMemory.get(),
                                                nonLinearSolverSens->get());
        break;
    case InternalSensitivityMethod::staggered1:
        status = CVodeSetNonlinearSolverSensStg1(solverMemory.get(),
                                                 nonLinearSolverSens->get());
        break;
    default:
        throw AmiException(
            "Unsupported internal sensitivity method selected: %d", ism);
    }

    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSolver::setNonLinearSolverSens");
}

void CVodeSolver::setNonLinearSolverB(const int which) const {
    int status = CVodeSetNonlinearSolverB(solverMemory.get(), which,
                                          nonLinearSolverB->get());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetNonlinearSolverB");
}

void CVodeSolver::setErrHandlerFn() const {
    int status =
        CVodeSetErrHandlerFn(solverMemory.get(), wrapErrHandlerFn,
                             reinterpret_cast<void*>(
                                 const_cast<CVodeSolver*>(this)));
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetErrHandlerFn");
}

void CVodeSolver::setUserData(Model *model) const {
    int status = CVodeSetUserData(solverMemory.get(), model);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetUserData");
}

void CVodeSolver::setUserDataB(const int which, Model *model) const {
    int status = CVodeSetUserDataB(solverMemory.get(), which, model);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetUserDataB");
}

void CVodeSolver::setMaxNumSteps(const long mxsteps) const {
    int status = CVodeSetMaxNumSteps(solverMemory.get(), mxsteps);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetMaxNumSteps");
}

void CVodeSolver::setStabLimDet(const int stldet) const {
    int status = CVodeSetStabLimDet(solverMemory.get(), stldet);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetStabLimDet");
}

void CVodeSolver::setStabLimDetB(const int which, const int stldet) const {
    int status = CVodeSetStabLimDetB(solverMemory.get(), which, stldet);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetStabLimDetB");
}

void CVodeSolver::setId(const Model * /*model*/) const {}

void CVodeSolver::setSuppressAlg(const bool /*flag*/) const {}

void CVodeSolver::resetState(void *ami_mem, const_N_Vector y0) const {

    auto cv_mem = static_cast<CVodeMem>(ami_mem);
    /* here we force the order in the next step to zero, and update the
     Nordsieck history array, this is largely copied from CVodeReInit with
     explanations from cvodes_impl.h
     */

    /* Set step parameters */

    /* current order */
    cv_mem->cv_q = 1;
    /* L = q + 1 */
    cv_mem->cv_L = 2;
    /* number of steps to wait before updating in q */
    cv_mem->cv_qwait = cv_mem->cv_L;
    /* last successful q value used                */
    cv_mem->cv_qu = 0;
    /* last successful h value used                */
    cv_mem->cv_hu = ZERO;
    /* tolerance scale factor                      */
    cv_mem->cv_tolsf = ONE;

    /* Initialize other integrator optional outputs */

    /* actual initial stepsize                     */
    cv_mem->cv_h0u = ZERO;
    /* step size to be used on the next step       */
    cv_mem->cv_next_h = ZERO;
    /* order to be used on the next step           */
    cv_mem->cv_next_q = 0;

    /* write updated state to Nordsieck history array  */
    N_VScale(ONE, y0, cv_mem->cv_zn[0]);
}

void CVodeSolver::reInitPostProcessF(const realtype tnext) const {
    reInitPostProcess(solverMemory.get(), &t, &x, tnext);
    forceReInitPostProcessF = false;
}

void CVodeSolver::reInitPostProcessB(const realtype tnext) const {
    realtype tBret;
    auto cv_mem = static_cast<CVodeMem>(solverMemory.get());
    auto ca_mem = cv_mem->cv_adj_mem;
    auto cvB_mem = ca_mem->cvB_mem;
    // loop over all backward problems
    while (cvB_mem != nullptr) {
        // store current backward problem in ca_mem to make it accessible in
        // adjoint rhs wrapper functions
        ca_mem->ca_bckpbCrt = cvB_mem;
        reInitPostProcess(static_cast<void *>(cvB_mem->cv_mem), &tBret, &xB,
                          tnext);
        cvB_mem->cv_tout = tBret;
        cvB_mem = cvB_mem->cv_next;
    }
    forceReInitPostProcessB = false;
}

void CVodeSolver::reInitPostProcess(void *ami_mem, realtype *t, AmiVector *yout,
                                    const realtype tout) const {
    auto cv_mem = static_cast<CVodeMem>(ami_mem);
    auto nst_tmp = cv_mem->cv_nst;
    cv_mem->cv_nst = 0;

    auto status = CVodeSetStopTime(cv_mem, tout);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetStopTime");

    status = CVode(ami_mem, tout, yout->getNVector(), t, CV_ONE_STEP);

    if (status != CV_SUCCESS)
        throw CvodeException(status, "reInitPostProcess");

    cv_mem->cv_nst = nst_tmp + 1;
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

void CVodeSolver::reInit(const realtype t0, const AmiVector &yy0,
                         const AmiVector & /*yp0*/) const {
    auto cv_mem = static_cast<CVodeMem>(solverMemory.get());
    cv_mem->cv_tn = t0;
    if (solverWasCalledF)
        forceReInitPostProcessF = true;
    x.copy(yy0);
    resetState(cv_mem, x.getNVector());
}

void CVodeSolver::sensReInit(const AmiVectorArray &yyS0,
                             const AmiVectorArray & /*ypS0*/) const {
    auto cv_mem = static_cast<CVodeMem>(solverMemory.get());
    /* Initialize znS[0] in the history array */
    for (int is = 0; is < nplist(); is++)
        cv_mem->cv_cvals[is] = ONE;
    if (solverWasCalledF)
        forceReInitPostProcessF = true;
    sx.copy(yyS0);
    int status = N_VScaleVectorArray(nplist(), cv_mem->cv_cvals,
                                     sx.getNVectorArray(), cv_mem->cv_znS[0]);
    if (status != CV_SUCCESS)
        throw CvodeException(CV_VECTOROP_ERR, "CVodeSensReInit");
}

void CVodeSolver::reInitB(const int which, const realtype tB0,
                          const AmiVector &yyB0,
                          const AmiVector & /*ypB0*/) const {
    auto cv_memB =
        static_cast<CVodeMem>(CVodeGetAdjCVodeBmem(solverMemory.get(), which));
    if (solverWasCalledB)
        forceReInitPostProcessB = true;
    cv_memB->cv_tn = tB0;
    xB.copy(yyB0);
    resetState(cv_memB, xB.getNVector());
}

void CVodeSolver::sensToggleOff() const {
    auto status = CVodeSensToggleOff(solverMemory.get());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSensToggleOff");
}

void CVodeSolver::quadReInitB(int which, const AmiVector &yQB0) const {
    auto cv_memB =
        static_cast<CVodeMem>(CVodeGetAdjCVodeBmem(solverMemory.get(), which));
    if (solverWasCalledB)
        forceReInitPostProcessB = true;
    xQB.copy(yQB0);
    N_VScale(ONE, xQB.getNVector(), cv_memB->cv_znQ[0]);
}

void CVodeSolver::setSensParams(const realtype *p, const realtype *pbar,
                                const int *plist) const {
    int status = CVodeSetSensParams(
        solverMemory.get(), const_cast<realtype *>(p),
        const_cast<realtype *>(pbar), const_cast<int *>(plist));
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetSensParams");
}

void CVodeSolver::getDky(realtype t, int k) const {
    int status = CVodeGetDky(solverMemory.get(), t, k, dky.getNVector());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetDky");
}

void CVodeSolver::getSens() const {
    realtype tDummy = 0;
    int status =
        CVodeGetSens(solverMemory.get(), &tDummy, sx.getNVectorArray());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetSens");
}

void CVodeSolver::getSensDky(const realtype t, const int k) const {
    int status =
        CVodeGetSensDky(solverMemory.get(), t, k, sx.getNVectorArray());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetSens");
}

void CVodeSolver::getDkyB(const realtype t, const int k,
                          const int which) const {
    int status = CVodeGetDky(CVodeGetAdjCVodeBmem(solverMemory.get(), which), t,
                             k, dky.getNVector());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetDkyB");
}

void CVodeSolver::getQuadB(int which) const {
    realtype tDummy = 0;
    int status =
        CVodeGetQuadB(solverMemory.get(), which, &tDummy, xQB.getNVector());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetQuadB");
}

void CVodeSolver::getQuad(const realtype t) const {
    int status =
    CVodeGetQuad(solverMemory.get(), &t, xQ.getNVector());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetQuad");
}

void CVodeSolver::getQuadDkyB(const realtype t, const int k, int which) const {
    int status =
        CVodeGetQuadDky(CVodeGetAdjCVodeBmem(solverMemory.get(), which), t, k,
                        xQB.getNVector());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetQuadDkyB");
}

void CVodeSolver::getQuadDky(const realtype t, const int k) const {
    int status =
    CVodeGetQuadDky(solverMemory.get(), t, k, xQ.getNVector());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetQuadDky");
}

void CVodeSolver::adjInit() const {
    int status;
    if (getAdjInitDone()) {
        status = CVodeAdjReInit(solverMemory.get());
    } else {
        status = CVodeAdjInit(solverMemory.get(), static_cast<int>(maxsteps),
                              static_cast<int>(interpType));
        setAdjInitDone();
    }
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeAdjInit");
}

void CVodeSolver::quadInit() const {
    int status;
    if (getQuadInitDone()) {
        status = CVodeQuadReInit(solverMemory.get(), xQ0.getNVector());
    } else {
        status = CVodeQuadInit(solverMemory.get(), fQdot_ss, xQ.getNVector());
        setQuadInitDone();
    }
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeQuadInit");
}

void CVodeSolver::allocateSolverB(int *which) const {
    if (!solverMemoryB.empty()) {
        *which = 0;
        return;
    }
    int status = CVodeCreateB(solverMemory.get(), static_cast<int>(lmm), which);
    if (*which + 1 > static_cast<int>(solverMemoryB.size()))
        solverMemoryB.resize(*which + 1);
    solverMemoryB.at(*which) =
        std::unique_ptr<void, std::function<void(void *)>>(
            getAdjBmem(solverMemory.get(), *which), [](void * /*ptr*/) {});
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeCreateB");
}

void CVodeSolver::setSStolerancesB(const int which, const realtype relTolB,
                                   const realtype absTolB) const {
    int status =
        CVodeSStolerancesB(solverMemory.get(), which, relTolB, absTolB);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSStolerancesB");
}

void CVodeSolver::quadSStolerancesB(const int which, const realtype reltolQB,
                                    const realtype abstolQB) const {
    int status =
        CVodeQuadSStolerancesB(solverMemory.get(), which, reltolQB, abstolQB);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeQuadSStolerancesB");
}

void CVodeSolver::quadSStolerances(const realtype reltolQB,
                                   const realtype abstolQB) const {
    int status =
    CVodeQuadSStolerances(solverMemory.get(),reltolQB, abstolQB);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeQuadSStolerances");
}

void CVodeSolver::getB(const int which) const {
    realtype tDummy = 0;
    int status = CVodeGetB(solverMemory.get(), which, &tDummy, xB.getNVector());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetB");
}

int CVodeSolver::solve(const realtype tout, const int itask) const {
    if (forceReInitPostProcessF)
        reInitPostProcessF(tout);
    int status = CVode(solverMemory.get(), tout, x.getNVector(), &t, itask);
    if (status < 0) // status > 0 is okay and is used for e.g. root return
        throw IntegrationFailure(status, t);
    solverWasCalledF = true;
    return status;
}

int CVodeSolver::solveF(const realtype tout, const int itask,
                        int *ncheckPtr) const {
    if (forceReInitPostProcessF)
        reInitPostProcessF(tout);
    int status =
        CVodeF(solverMemory.get(), tout, x.getNVector(), &t, itask, ncheckPtr);
    if (status < 0) // status > 0 is okay and is used for e.g. root return
        throw IntegrationFailure(status, t);
    solverWasCalledF = true;
    return status;
}

void CVodeSolver::solveB(const realtype tBout, const int itaskB) const {
    if (forceReInitPostProcessB)
        reInitPostProcessB(tBout);
    int status = CVodeB(solverMemory.get(), tBout, itaskB);
    if (status != CV_SUCCESS)
        throw IntegrationFailureB(status, tBout);
    solverWasCalledB = true;
}

void CVodeSolver::setMaxNumStepsB(const int which, const long mxstepsB) const {
    int status = CVodeSetMaxNumStepsB(solverMemory.get(), which, mxstepsB);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetMaxNumStepsB");
}

void CVodeSolver::diag() const {
    int status = CVDiag(solverMemory.get());
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVDiag");
}

void CVodeSolver::diagB(const int which) const {
    int status = CVDiagB(solverMemory.get(), which);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVDiagB");
}

void CVodeSolver::getNumSteps(const void *ami_mem, long int *numsteps) const {
    int status = CVodeGetNumSteps(const_cast<void *>(ami_mem), numsteps);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetNumSteps");
}

void CVodeSolver::getNumRhsEvals(const void *ami_mem,
                                 long int *numrhsevals) const {
    int status = CVodeGetNumRhsEvals(const_cast<void *>(ami_mem), numrhsevals);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetNumRhsEvals");
}

void CVodeSolver::getNumErrTestFails(const void *ami_mem,
                                     long int *numerrtestfails) const {
    int status =
        CVodeGetNumErrTestFails(const_cast<void *>(ami_mem), numerrtestfails);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetNumErrTestFails");
}

void CVodeSolver::getNumNonlinSolvConvFails(
    const void *ami_mem, long int *numnonlinsolvconvfails) const {
    int status = CVodeGetNumNonlinSolvConvFails(const_cast<void *>(ami_mem),
                                                numnonlinsolvconvfails);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetNumNonlinSolvConvFails");
}

void CVodeSolver::getLastOrder(const void *ami_mem, int *order) const {
    int status = CVodeGetLastOrder(const_cast<void *>(ami_mem), order);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeGetLastOrder");
}

void *CVodeSolver::getAdjBmem(void *ami_mem, int which) const {
    return CVodeGetAdjCVodeBmem(ami_mem, which);
}

void CVodeSolver::calcIC(const realtype /*tout1*/) const {};

void CVodeSolver::calcICB(const int /*which*/, const realtype /*tout1*/) const {};

void CVodeSolver::setStopTime(const realtype tstop) const {
    int status = CVodeSetStopTime(solverMemory.get(), tstop);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeSetStopTime");
}


void CVodeSolver::turnOffRootFinding() const {
    int status = CVodeRootInit(solverMemory.get(), 0, nullptr);
    if (status != CV_SUCCESS)
        throw CvodeException(status, "CVodeRootInit");
}


const Model *CVodeSolver::getModel() const {
    if (!solverMemory)
        throw AmiException("Solver has not been allocated, information is not "
                           "available");
    auto cv_mem = static_cast<CVodeMem>(solverMemory.get());
    return static_cast<Model *>(cv_mem->cv_user_data);
}


/**
 * @brief Jacobian of xdot with respect to states x
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
static int fJ(realtype t, N_Vector x, N_Vector xdot, SUNMatrix J,
                    void *user_data, N_Vector /*tmp1*/, N_Vector /*tmp2*/,
                    N_Vector /*tmp3*/) {
    auto model = static_cast<Model_ODE *>(user_data);
    model->fJ(t, x, xdot, J);
    return model->checkFinite(gsl::make_span(J), "Jacobian");
}


/**
 * @brief Jacobian of xBdot with respect to adjoint state xB
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
static int fJB(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot,
                     SUNMatrix JB, void *user_data, N_Vector /*tmp1B*/,
                     N_Vector /*tmp2B*/, N_Vector /*tmp3B*/) {
    auto model = static_cast<Model_ODE *>(user_data);
    model->fJB(t, x, xB, xBdot, JB);
    return model->checkFinite(gsl::make_span(JB), "Jacobian");
}


/**
 * @brief J in sparse form (for sparse solvers from the SuiteSparse Package)
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
static int fJSparse(realtype t, N_Vector x, N_Vector /*xdot*/,
                          SUNMatrix J, void *user_data, N_Vector /*tmp1*/,
                          N_Vector /*tmp2*/, N_Vector /*tmp3*/) {
    auto model = static_cast<Model_ODE *>(user_data);
    model->fJSparse(t, x, J);
    return model->checkFinite(gsl::make_span(J), "Jacobian");
}


/**
 * @brief JB in sparse form (for sparse solvers from the SuiteSparse Package)
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
static int fJSparseB(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot,
                           SUNMatrix JB, void *user_data, N_Vector /*tmp1B*/,
                           N_Vector /*tmp2B*/, N_Vector /*tmp3B*/) {
    auto model = static_cast<Model_ODE *>(user_data);
    model->fJSparseB(t, x, xB, xBdot, JB);
    return model->checkFinite(gsl::make_span(JB), "Jacobian");
}


/**
 * @brief J in banded form (for banded solvers)
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
static int fJBand(realtype t, N_Vector x, N_Vector xdot, SUNMatrix J,
           void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    return fJ(t, x, xdot, J, user_data, tmp1, tmp2, tmp3);
}


/**
 * @brief JB in banded form (for banded solvers)
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
static int fJBandB(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot,
                         SUNMatrix JB, void *user_data, N_Vector tmp1B,
                         N_Vector tmp2B, N_Vector tmp3B) {
    return fJB(t, x, xB, xBdot, JB, user_data, tmp1B, tmp2B, tmp3B);
}


/**
 * @brief Matrix vector product of J with a vector v (for iterative solvers)
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
static int fJv(N_Vector v, N_Vector Jv, realtype t, N_Vector x,
        N_Vector /*xdot*/, void *user_data, N_Vector /*tmp*/) {
    auto model = static_cast<Model_ODE *>(user_data);
    model->fJv(v, Jv, t, x);
    return model->checkFinite(gsl::make_span(Jv), "Jacobian");
}


/**
 * @brief Matrix vector product of JB with a vector v (for iterative solvers)
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
static int fJvB(N_Vector vB, N_Vector JvB, realtype t, N_Vector x,
         N_Vector xB, N_Vector /*xBdot*/, void *user_data,
         N_Vector /*tmpB*/) {
    auto model = static_cast<Model_ODE *>(user_data);
    model->fJvB(vB, JvB, t, x, xB);
    return model->checkFinite(gsl::make_span(JvB), "Jacobian");
}


/**
 * @brief Event trigger function for events
 * @param t timepoint
 * @param x Vector with the states
 * @param root array with root function values
 * @param user_data object with user input @type Model_ODE
 * @return status flag indicating successful execution
 */
static int froot(realtype t, N_Vector x, realtype *root,
                       void *user_data) {
    auto model = static_cast<Model_ODE *>(user_data);
    model->froot(t, x, gsl::make_span<realtype>(root, model->ne));
    return model->checkFinite(gsl::make_span<realtype>(root, model->ne),
                              "root function");
}


/**
 * @brief residual function of the ODE
 * @param t timepoint
 * @param x Vector with the states
 * @param xdot Vector with the right hand side
 * @param user_data object with user input @type Model_ODE
 * @return status flag indicating successful execution
 */
static int fxdot(realtype t, N_Vector x, N_Vector xdot, void *user_data) {
    auto model = static_cast<Model_ODE *>(user_data);

    if (t > 1e200 && !model->checkFinite(gsl::make_span(x), "fxdot")) {
        /* when t is large (typically ~1e300), CVODES may pass all NaN x
           to fxdot from which we typically cannot recover. To save time
           on normal execution, we do not always want to check finiteness
           of x, but only do so when t is large and we expect problems. */
        return AMICI_UNRECOVERABLE_ERROR;
    }

    model->fxdot(t, x, xdot);
    return model->checkFinite(gsl::make_span(xdot), "fxdot");
}


/**
 * @brief Right hand side of differential equation for adjoint state xB
 * @param t timepoint
 * @param x Vector with the states
 * @param xB Vector with the adjoint states
 * @param xBdot Vector with the adjoint right hand side
 * @param user_data object with user input @type Model_ODE
 * @return status flag indicating successful execution
 */
static int fxBdot(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot,
                        void *user_data) {
    auto model = static_cast<Model_ODE *>(user_data);
    model->fxBdot(t, x, xB, xBdot);
    return model->checkFinite(gsl::make_span(xBdot), "fxBdot");
}


/**
 * @brief Right hand side of integral equation for quadrature states qB
 * @param t timepoint
 * @param x Vector with the states
 * @param xB Vector with the adjoint states
 * @param qBdot Vector with the adjoint quadrature right hand side
 * @param user_data pointer to temp data object
 * @return status flag indicating successful execution
 */
static int fqBdot(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot,
                        void *user_data) {
    auto model = static_cast<Model_ODE *>(user_data);
    model->fqBdot(t, x, xB, qBdot);
    return model->checkFinite(gsl::make_span(qBdot), "qBdot");
}


/**
 * @brief Right hand side of differential equation for state sensitivities sx
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
static int fsxdot(int /*Ns*/, realtype t, N_Vector x, N_Vector /*xdot*/,
                        int ip, N_Vector sx, N_Vector sxdot, void *user_data,
                        N_Vector /*tmp1*/, N_Vector /*tmp2*/) {
    auto model = static_cast<Model_ODE *>(user_data);
    model->fsxdot(t, x, ip, sx, sxdot);
    return model->checkFinite(gsl::make_span(sxdot), "sxdot");
}


bool operator==(const CVodeSolver &a, const CVodeSolver &b) {
    return static_cast<Solver const &>(a) == static_cast<Solver const &>(b);
}

} // namespace amici
