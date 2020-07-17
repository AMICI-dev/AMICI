#include "amici/solver_idas.h"

#include "amici/exception.h"
#include "amici/misc.h"
#include "amici/model_dae.h"
#include "amici/sundials_linsol_wrapper.h"

#include <idas/idas.h>
#include <idas/idas_impl.h>

#include <amd.h>
#include <btf.h>
#include <colamd.h>
#include <klu.h>

#define ONE RCONST(1.0)

namespace amici {

/*
 * The following static members are callback function to CVODES.
 * Their signatures must not be changes.
 */

static int fxdot(realtype t, N_Vector x, N_Vector dx, N_Vector xdot,
                 void *user_data);

static int fJ(realtype t, realtype cj, N_Vector x, N_Vector dx,
              N_Vector xdot, SUNMatrix J, void *user_data, N_Vector tmp1,
              N_Vector tmp2, N_Vector tmp3);

static int fJSparse(realtype t, realtype cj, N_Vector x, N_Vector dx,
                    N_Vector xdot, SUNMatrix J, void *user_data,
                    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int fJB(realtype t, realtype cj, N_Vector x, N_Vector dx,
               N_Vector xB, N_Vector dxB, N_Vector xBdot, SUNMatrix JB,
               void *user_data, N_Vector tmp1B, N_Vector tmp2B,
               N_Vector tmp3B);

static int fJSparseB(realtype t, realtype cj, N_Vector x, N_Vector dx,
                     N_Vector xB, N_Vector dxB, N_Vector xBdot,
                     SUNMatrix JB, void *user_data, N_Vector tmp1B,
                     N_Vector tmp2B, N_Vector tmp3B);

static int fJBand(realtype t, realtype cj, N_Vector x, N_Vector dx,
                  N_Vector xdot, SUNMatrix J, void *user_data,
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int fJBandB(realtype t, realtype cj, N_Vector x, N_Vector dx,
                   N_Vector xB, N_Vector dxB, N_Vector xBdot, SUNMatrix JB,
                   void *user_data, N_Vector tmp1B, N_Vector tmp2B,
                   N_Vector tmp3B);

static int fJv(realtype t, N_Vector x, N_Vector dx, N_Vector xdot,
               N_Vector v, N_Vector Jv, realtype cj, void *user_data,
               N_Vector tmp1, N_Vector tmp2);

static int fJvB(realtype t, N_Vector x, N_Vector dx, N_Vector xB,
                N_Vector dxB, N_Vector xBdot, N_Vector vB, N_Vector JvB,
                realtype cj, void *user_data, N_Vector tmpB1,
                N_Vector tmpB2);

static int froot(realtype t, N_Vector x, N_Vector dx, realtype *root,
                 void *user_data);

static int fxBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB,
                  N_Vector dxB, N_Vector xBdot, void *user_data);

static int fqBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB,
                  N_Vector dxB, N_Vector qBdot, void *user_data);

static int fxBdot_ss(realtype t, N_Vector xB, N_Vector dxB, N_Vector xBdot,
                     void *user_data);

static int fqBdot_ss(realtype t, N_Vector xB, N_Vector dxB, N_Vector qBdot,
                     void *user_data);

static int fJSparseB_ss(realtype t, realtype cj, N_Vector x, N_Vector dx,
                        N_Vector xBdot, SUNMatrix JB, void *user_data,
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int fsxdot(int Ns, realtype t, N_Vector x, N_Vector dx,
                  N_Vector xdot, N_Vector *sx, N_Vector *sdx,
                  N_Vector *sxdot, void *user_data, N_Vector tmp1,
                  N_Vector tmp2, N_Vector tmp3);


/* Function implementations */

void IDASolver::init(const realtype t0, const AmiVector &x0,
                     const AmiVector &dx0) const {
    int status;
    solverWasCalledF = false;
    t = t0;
    x = x0;
    dx = dx0;
    if (getInitDone()) {
        status =
            IDAReInit(solverMemory.get(), t, x.getNVector(), dx.getNVector());
    } else {
        status = IDAInit(solverMemory.get(), fxdot, t, x.getNVector(),
                         dx.getNVector());
        setInitDone();
    }
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDAInit");
}

void IDASolver::initSteadystate(const realtype t0, const AmiVector &x0,
                                const AmiVector &dx0) const {
    /* We need to set the steadystate rhs function. SUndials doesn't have this
       in its public api, so we have to change it in the solver memory,
       as re-calling init would unset solver settings. */
    auto ida_mem = static_cast<IDAMem>(solverMemory.get());
    ida_mem->ida_res = fxBdot_ss;
}

void IDASolver::sensInit1(const AmiVectorArray &sx0,
                          const AmiVectorArray &sdx0) const {
    int status = IDA_SUCCESS;
    sx = sx0;
    sdx = sdx0;
    if (getSensitivityMethod() == SensitivityMethod::forward && nplist() > 0) {
        if (getSensInitDone()) {
            status =
                IDASensReInit(solverMemory.get(),
                              static_cast<int>(getInternalSensitivityMethod()),
                              sx.getNVectorArray(), sdx.getNVectorArray());
        } else {
            status = IDASensInit(
                solverMemory.get(), nplist(),
                static_cast<int>(getInternalSensitivityMethod()), fsxdot,
                sx.getNVectorArray(), sdx.getNVectorArray());
            setSensInitDone();
        }
    }
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDASensInit");
}

void IDASolver::binit(const int which, const realtype tf, const AmiVector &xB0,
                      const AmiVector &dxB0) const {
    int status;
    xB = xB0;
    dxB = dxB0;
    if (getInitDoneB(which))
        status = IDAReInitB(solverMemory.get(), which, tf, xB.getNVector(),
                            dxB.getNVector());
    else {

        status = IDAInitB(solverMemory.get(), which, fxBdot, tf,
                          xB.getNVector(), dxB.getNVector());
        setInitDoneB(which);
    }
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDAInitB");
}

void IDASolver::qbinit(const int which, const AmiVector &xQB0) const {
    int status;
    xQB.copy(xQB0);
    if (getQuadInitDoneB(which))
        status = IDAQuadReInitB(solverMemory.get(), which, xQB.getNVector());
    else {
        status =
            IDAQuadInitB(solverMemory.get(), which, fqBdot, xQB.getNVector());
        setQuadInitDoneB(which);
    }
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDAQuadInitB");
}

void IDASolver::rootInit(int ne) const {
    int status = IDARootInit(solverMemory.get(), ne, froot);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDARootInit");
}

void IDASolver::setDenseJacFn() const {
    int status = IDASetJacFn(solverMemory.get(), fJ);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDADlsSetDenseJacFn");
}

void IDASolver::setSparseJacFn() const {
    int status = IDASetJacFn(solverMemory.get(), fJSparse);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDASlsSetSparseJacFn");
}

void IDASolver::setBandJacFn() const {
    int status = IDASetJacFn(solverMemory.get(), fJBand);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDADlsSetBandJacFn");
}

void IDASolver::setJacTimesVecFn() const {
    int status = IDASetJacTimes(solverMemory.get(), nullptr, fJv);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDASpilsSetJacTimesVecFn");
}

void IDASolver::setDenseJacFnB(const int which) const {
    int status = IDASetJacFnB(solverMemory.get(), which, fJB);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDADlsSetDenseJacFnB");
}

void IDASolver::setSparseJacFnB(const int which) const {
    int status = IDASetJacFnB(solverMemory.get(), which, fJSparseB);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDASlsSetSparseJacFnB");
}

void IDASolver::setBandJacFnB(const int which) const {
    int status = IDASetJacFnB(solverMemory.get(), which, fJBandB);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDADlsSetBandJacFnB");
}

void IDASolver::setJacTimesVecFnB(const int which) const {
    int status = IDASetJacTimesB(solverMemory.get(), which, nullptr, fJvB);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDASpilsSetJacTimesVecFnB");
}

void IDASolver::setSparseJacFn_ss() const {
    int status = IDASetJacFn(solverMemory.get(), fJSparseB_ss);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDASetJacFn");
}

Solver *IDASolver::clone() const { return new IDASolver(*this); }

void IDASolver::allocateSolver() const {
    if (!solverMemory)
        solverMemory = std::unique_ptr<void, std::function<void(void *)>>(
            IDACreate(), [](void *ptr) { IDAFree(&ptr); });
}

void IDASolver::setSStolerances(const realtype rtol,
                                const realtype atol) const {
    int status = IDASStolerances(solverMemory.get(), rtol, atol);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDASStolerances");
}
void IDASolver::setSensSStolerances(const realtype rtol,
                                    const realtype *atol) const {
    int status = IDASensSStolerances(solverMemory.get(), rtol,
                                     const_cast<realtype *>(atol));
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDASensEEtolerances");
}

void IDASolver::setSensErrCon(const bool error_corr) const {
    int status = IDASetSensErrCon(solverMemory.get(), error_corr);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDASetSensErrCon");
}

void IDASolver::setQuadErrConB(const int which, const bool flag) const {
    int status = IDASetQuadErrConB(solverMemory.get(), which, flag);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDASetQuadErrConB");
}

void IDASolver::setQuadErrCon(const bool flag) const {
    int status = IDASetQuadErrCon(solverMemory.get(), flag);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDASetQuadErrCon");
}

void IDASolver::getRootInfo(int *rootsfound) const {
    int status = IDAGetRootInfo(solverMemory.get(), rootsfound);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDAGetRootInfo");
}

void IDASolver::setErrHandlerFn() const {
    int status =
        IDASetErrHandlerFn(solverMemory.get(), wrapErrHandlerFn,
                           reinterpret_cast<void*>(
                               const_cast<IDASolver*>(this)));
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDASetErrHandlerFn");
}

void IDASolver::setUserData(Model *model) const {
    int status = IDASetUserData(solverMemory.get(), model);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDASetUserData");
}

void IDASolver::setUserDataB(int which, Model *model) const {
    int status = IDASetUserDataB(solverMemory.get(), which, model);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDASetUserDataB");
}

void IDASolver::setMaxNumSteps(const long int mxsteps) const {
    int status = IDASetMaxNumSteps(solverMemory.get(), mxsteps);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDASetMaxNumSteps");
}

void IDASolver::setStabLimDet(const int /*stldet*/) const {}

void IDASolver::setStabLimDetB(const int /*which*/,
                               const int /*stldet*/) const {}

void IDASolver::setId(const Model *model) const {

    N_Vector id = N_VMake_Serial(model->nx_solver,
                                 const_cast<realtype *>(model->idlist.data()));

    int status = IDASetId(solverMemory.get(), id);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDASetMaxNumSteps");

    N_VDestroy_Serial(id);
}

void IDASolver::setSuppressAlg(const bool flag) const {
    int status = IDASetSuppressAlg(solverMemory.get(), flag);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDASetSuppressAlg");
}

void IDASolver::resetState(void *ami_mem, const_N_Vector yy0,
                           const_N_Vector yp0) const {

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
}

void IDASolver::reInitPostProcessF(const realtype tnext) const {
    reInitPostProcess(solverMemory.get(), &t, &x, &dx, tnext);
}

void IDASolver::reInitPostProcessB(const realtype tnext) const {
    realtype tBret;
    auto ida_mem = static_cast<IDAMem>(solverMemory.get());
    auto idaadj_mem = ida_mem->ida_adj_mem;
    auto idaB_mem = idaadj_mem->IDAB_mem;
    // loop over all backward problems
    while (idaB_mem != nullptr) {
        // store current backward problem in ca_mem to make it accessible in
        // adjoint rhs wrapper functions
        idaadj_mem->ia_bckpbCrt = idaB_mem;
        reInitPostProcess(static_cast<void *>(idaB_mem->IDA_mem), &tBret, &xB,
                          &dxB, tnext);
        // idaB_mem->ida_tout = tBret;
        idaB_mem = idaB_mem->ida_next;
    }
    forceReInitPostProcessB = false;
}

void IDASolver::reInitPostProcess(void *ami_mem, realtype *t,
                                  AmiVector *yout, AmiVector *ypout,
                                  realtype tout) const {
    auto ida_mem = static_cast<IDAMem>(ami_mem);
    auto nst_tmp = ida_mem->ida_nst;
    ida_mem->ida_nst = 0;

    auto status = IDASetStopTime(ida_mem, tout);
    if(status != IDA_SUCCESS)
        throw IDAException(status, "CVodeSetStopTime");

    status = IDASolve(ami_mem, tout, t, yout->getNVector(), ypout->getNVector(),
                      IDA_ONE_STEP);

    if(status != IDA_SUCCESS)
        throw IDAException(status, "reInitPostProcess");

    ida_mem->ida_nst = nst_tmp+1;
    if (ida_mem->ida_adjMallocDone == SUNTRUE) {
        /* add new step to history array, this is copied from CVodeF */
        auto ia_mem = ida_mem->ida_adj_mem;
        auto dt_mem = ia_mem->dt_mem;

        if (ida_mem->ida_nst % ia_mem->ia_nsteps == 0) {
            /* currently not implemented, we should never get here as we
             limit cv_mem->cv_nst < ca_mem->ca_nsteps, keeping this for
             future regression */
            throw IDAException(AMICI_ERROR, "reInitPostProcess");
        }

        /* Load next point in dt_mem */
        dt_mem[ida_mem->ida_nst % ia_mem->ia_nsteps]->t = *t;
        ia_mem->ia_storePnt(ida_mem,
                            dt_mem[ida_mem->ida_nst % ia_mem->ia_nsteps]);

        /* Set t1 field of the current ckeck point structure
         for the case in which there will be no future
         check points */
        ia_mem->ck_mem->ck_t1 = *t;

        /* tfinal is now set to *tret */
        ia_mem->ia_tfinal = *t;
    }
}

void IDASolver::reInit(const realtype t0, const AmiVector &yy0,
                       const AmiVector &yp0) const {

    auto ida_mem = static_cast<IDAMem>(solverMemory.get());
    ida_mem->ida_tn = t0;
    if (solverWasCalledF)
        forceReInitPostProcessF = true;
    x.copy(yy0);
    dx.copy(yp0);
    resetState(ida_mem, x.getNVector(), xB.getNVector());
}

void IDASolver::sensReInit(const AmiVectorArray &yyS0,
                           const AmiVectorArray &ypS0) const {
    auto ida_mem = static_cast<IDAMem>(solverMemory.get());
    /* Initialize znS[0] in the history array */
    for (int is = 0; is < nplist(); is++)
        ida_mem->ida_cvals[is] = ONE;
    if (solverWasCalledF)
        forceReInitPostProcessF = true;
    sx.copy(yyS0);
    sdx.copy(ypS0);
    auto status =
        N_VScaleVectorArray(nplist(), ida_mem->ida_cvals, sx.getNVectorArray(),
                            ida_mem->ida_phiS[0]);
    if (status != IDA_SUCCESS)
        throw IDAException(IDA_VECTOROP_ERR, "IDASensReInit");
    status = N_VScaleVectorArray(nplist(), ida_mem->ida_cvals,
                                 sdx.getNVectorArray(), ida_mem->ida_phiS[1]);
    if (status != IDA_SUCCESS)
        throw IDAException(IDA_VECTOROP_ERR, "IDASensReInit");
}

void IDASolver::sensToggleOff() const {
    auto status = IDASensToggleOff(solverMemory.get());
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDASensToggleOff");
}

void IDASolver::reInitB(const int which, const realtype tB0,
                        const AmiVector &yyB0, const AmiVector &ypB0) const {

    auto ida_memB =
        static_cast<IDAMem>(IDAGetAdjIDABmem(solverMemory.get(), which));
    if (solverWasCalledB)
        forceReInitPostProcessB = true;
    ida_memB->ida_tn = tB0;
    xB.copy(yyB0);
    dxB.copy(ypB0);
    resetState(ida_memB, xB.getNVector(), dxB.getNVector());
}

void IDASolver::quadReInitB(const int which, const AmiVector &yQB0) const {
    auto ida_memB =
        static_cast<IDAMem>(IDAGetAdjIDABmem(solverMemory.get(), which));
    if (solverWasCalledB)
        forceReInitPostProcessB = true;
    xQB.copy(yQB0);
    N_VScale(ONE, xQB.getNVector(), ida_memB->ida_phiQ[0]);
}

void IDASolver::setSensParams(const realtype *p, const realtype *pbar,
                              const int *plist) const {
    int status = IDASetSensParams(solverMemory.get(), const_cast<realtype *>(p),
                                  const_cast<realtype *>(pbar),
                                  const_cast<int *>(plist));
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDASetSensParams");
}

void IDASolver::getDky(const realtype t, const int k) const {
    int status = IDAGetDky(solverMemory.get(), t, k, dky.getNVector());
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDAGetDky");
}

void IDASolver::getSens() const {
    realtype tDummy = 0;
    int status = IDAGetSens(solverMemory.get(), &tDummy, sx.getNVectorArray());
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDAGetSens");
}

void IDASolver::getSensDky(const realtype t, const int k) const {
    int status = IDAGetSensDky(solverMemory.get(), t, k, sx.getNVectorArray());
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDAGetSens");
}

void IDASolver::getB(const int which) const {
    realtype tDummy = 0;
    int status = IDAGetB(solverMemory.get(), which, &tDummy, xB.getNVector(),
                         dxB.getNVector());
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDAGetB");
}

void IDASolver::getDkyB(const realtype t, int k, const int which) const {
    int status = IDAGetDky(IDAGetAdjIDABmem(solverMemory.get(), which), t, k,
                           dky.getNVector());
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDAGetB");
}
void IDASolver::getQuadB(int which) const {
    realtype tDummy = 0;
    int status =
        IDAGetQuadB(solverMemory.get(), which, &tDummy, xQB.getNVector());
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDAGetQuadB");
}

void IDASolver::getQuad(realtype &t) const {
    int status = IDAGetQuad(solverMemory.get(), &t, xQ.getNVector());
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDAGetQuad");
}

void IDASolver::getQuadDkyB(const realtype t, int k, const int which) const {
    int status = IDAGetQuadDky(IDAGetAdjIDABmem(solverMemory.get(), which), t,
                               k, xQB.getNVector());
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDAGetB");
}

void IDASolver::getQuadDky(const realtype t, const int k) const {
    int status = IDAGetQuadDky(solverMemory.get(), t, k, xQ.getNVector());
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDAGetQuadDky");
}

void IDASolver::adjInit() const {
    int status;
    if (getAdjInitDone()) {
        status = IDAAdjReInit(solverMemory.get());
    } else {
        status = IDAAdjInit(solverMemory.get(), static_cast<int>(maxsteps),
                            static_cast<int>(interpType));
        setAdjInitDone();
    }
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDAAdjInit");
}

void IDASolver::quadInit(const AmiVector &xQ0) const {
    int status;
    xQ.copy(xQ0);
    if (getQuadInitDone()) {
        status = IDAQuadReInit(solverMemory.get(), xQ0.getNVector());
    } else {
        status = IDAQuadInit(solverMemory.get(), fqBdot_ss, xQ.getNVector());
        setQuadInitDone();
    }
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDAQuadInit");
}

void IDASolver::allocateSolverB(int *which) const {
    if (!solverMemoryB.empty()) {
        *which = 0;
        return;
    }
    int status = IDACreateB(solverMemory.get(), which);
    if (*which + 1 > static_cast<int>(solverMemoryB.size()))
        solverMemoryB.resize(*which + 1);
    solverMemoryB.at(*which) =
        std::unique_ptr<void, std::function<void(void *)>>(
            getAdjBmem(solverMemory.get(), *which), [](void * /*ptr*/) {});
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDACreateB");
}

void IDASolver::setSStolerancesB(const int which, const realtype relTolB,
                                 const realtype absTolB) const {
    int status = IDASStolerancesB(solverMemory.get(), which, relTolB, absTolB);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDASStolerancesB");
}

void IDASolver::quadSStolerancesB(const int which, const realtype reltolQB,
                                  const realtype abstolQB) const {
    int status =
        IDAQuadSStolerancesB(solverMemory.get(), which, reltolQB, abstolQB);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDAQuadSStolerancesB");
}

void IDASolver::quadSStolerances(const realtype reltolQB,
                                 const realtype abstolQB) const {
    int status = IDAQuadSStolerances(solverMemory.get(), reltolQB, abstolQB);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDAQuadSStolerances");
}


int IDASolver::solve(const realtype tout, const int itask) const {
    if (forceReInitPostProcessF)
        reInitPostProcessF(tout);
    int status = IDASolve(solverMemory.get(), tout, &t, x.getNVector(),
                          dx.getNVector(), itask);
    solverWasCalledF = true;
    if (status < 0) // status > 0 is okay and is used for e.g. root return
        throw IntegrationFailure(status, t);
    return status;
}

int IDASolver::solveF(const realtype tout, const int itask,
                      int *ncheckPtr) const {
    if (forceReInitPostProcessF)
        reInitPostProcessF(tout);
    int status = IDASolveF(solverMemory.get(), tout, &t, x.getNVector(),
                           xB.getNVector(), itask, ncheckPtr);
    solverWasCalledF = true;
    if (status < 0) // status > 0 is okay and is used for e.g. root return
        throw IntegrationFailure(status, t);
    return status;
}

void IDASolver::solveB(const realtype tBout, const int itaskB) const {
    if (forceReInitPostProcessB)
        reInitPostProcessB(tBout);
    int status = IDASolveB(solverMemory.get(), tBout, itaskB);
    solverWasCalledB = true;
    if (status != IDA_SUCCESS)
        throw IntegrationFailure(status, tBout);
}

void IDASolver::setMaxNumStepsB(const int which,
                                const long int mxstepsB) const {
    int status = IDASetMaxNumStepsB(solverMemory.get(), which, mxstepsB);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDASetMaxNumStepsB");
}

void IDASolver::diag() const {
    throw AmiException("Diag Solver was not implemented for DAEs");
}

void IDASolver::diagB(const int /*which*/) const {
    throw AmiException("Diag Solver was not implemented for DAEs");
}

void IDASolver::getNumSteps(const void *ami_mem, long int *numsteps) const {
    int status = IDAGetNumSteps(const_cast<void *>(ami_mem), numsteps);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDAGetNumSteps");
}

void IDASolver::getNumRhsEvals(const void *ami_mem,
                               long int *numrhsevals) const {
    int status = IDAGetNumResEvals(const_cast<void *>(ami_mem), numrhsevals);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDAGetNumResEvals");
}

void IDASolver::getNumErrTestFails(const void *ami_mem,
                                   long int *numerrtestfails) const {
    int status =
        IDAGetNumErrTestFails(const_cast<void *>(ami_mem), numerrtestfails);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDAGetNumErrTestFails");
}

void IDASolver::getNumNonlinSolvConvFails(
    const void *ami_mem, long int *numnonlinsolvconvfails) const {
    int status = IDAGetNumNonlinSolvConvFails(const_cast<void *>(ami_mem),
                                              numnonlinsolvconvfails);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDAGetNumNonlinSolvConvFails");
}

void IDASolver::getLastOrder(const void *ami_mem, int *order) const {
    int status = IDAGetLastOrder(const_cast<void *>(ami_mem), order);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDAGetLastOrder");
}

void *IDASolver::getAdjBmem(void *ami_mem, int which) const {
    return IDAGetAdjIDABmem(ami_mem, which);
}

void IDASolver::calcIC(realtype tout1) const {
    int status = IDACalcIC(solverMemory.get(), IDA_YA_YDP_INIT, tout1);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDACalcIC");
    status =
        IDAGetConsistentIC(solverMemory.get(), x.getNVector(), dx.getNVector());
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDACalcIC");
}

void IDASolver::calcICB(const int which, const realtype tout1) const {
    int status = IDACalcICB(solverMemory.get(), which, tout1, xB.getNVector(),
                            dxB.getNVector());
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDACalcICB");
}

void IDASolver::setStopTime(const realtype tstop) const {
    int status = IDASetStopTime(solverMemory.get(), tstop);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDASetStopTime");
}

void IDASolver::turnOffRootFinding() const {
    int status = IDARootInit(solverMemory.get(), 0, nullptr);
    if (status != IDA_SUCCESS)
        throw IDAException(status, "IDARootInit");
}

const Model *IDASolver::getModel() const {
    if (!solverMemory)
        throw AmiException(
            "Solver has not been allocated, information is not available");
    auto ida_mem = static_cast<IDAMem>(solverMemory.get());
    return static_cast<Model *>(ida_mem->ida_user_data);
}

void IDASolver::setLinearSolver() const {
    int status = IDASetLinearSolver(solverMemory.get(), linearSolver->get(),
                                    linearSolver->getMatrix());
    if (status != IDA_SUCCESS)
        throw IDAException(status, "setLinearSolver");
}

void IDASolver::setLinearSolverB(const int which) const {
    int status =
        IDASetLinearSolverB(solverMemoryB[which].get(), which,
                            linearSolverB->get(), linearSolverB->getMatrix());
    if (status != IDA_SUCCESS)
        throw IDAException(status, "setLinearSolverB");
}

void IDASolver::setNonLinearSolver() const {
    int status =
        IDASetNonlinearSolver(solverMemory.get(), nonLinearSolver->get());
    if (status != IDA_SUCCESS)
        throw CvodeException(status, "CVodeSetNonlinearSolver");
}

void IDASolver::setNonLinearSolverSens() const {
    if (getSensitivityOrder() < SensitivityOrder::first)
        return;
    if (getSensitivityMethod() != SensitivityMethod::forward)
        return;

    int status = IDA_SUCCESS;

    switch (ism) {
    case InternalSensitivityMethod::staggered:
        status = IDASetNonlinearSolverSensStg(solverMemory.get(),
                                              nonLinearSolverSens->get());
        break;
    case InternalSensitivityMethod::simultaneous:
        status = IDASetNonlinearSolverSensSim(solverMemory.get(),
                                              nonLinearSolverSens->get());
        break;
    case InternalSensitivityMethod::staggered1:
    default:
        throw AmiException(
            "Unsupported internal sensitivity method selected: %d", ism);
    }

    if (status != IDA_SUCCESS)
        throw CvodeException(status, "CVodeSolver::setNonLinearSolverSens");
}

void IDASolver::setNonLinearSolverB(int which) const {
    int status = IDASetNonlinearSolverB(solverMemory.get(), which,
                                        nonLinearSolverB->get());
    if (status != IDA_SUCCESS)
        throw CvodeException(status, "CVodeSetNonlinearSolverB");
}

/**
 * @brief Jacobian of xdot with respect to states x
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
int fJ(realtype t, realtype cj, N_Vector x, N_Vector dx,
                  N_Vector xdot, SUNMatrix J, void *user_data,
                  N_Vector /*tmp1*/, N_Vector /*tmp2*/, N_Vector /*tmp3*/) {

    auto model = static_cast<Model_DAE *>(user_data);
    model->fJ(t, cj, x, dx, xdot, J);
    return model->checkFinite(gsl::make_span(J), "Jacobian");
}

/**
 * @brief Jacobian of xBdot with respect to adjoint state xB
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
int fJB(realtype t, realtype cj, N_Vector x, N_Vector dx,
                   N_Vector xB, N_Vector dxB, N_Vector /*xBdot*/, SUNMatrix JB,
                   void *user_data, N_Vector /*tmp1B*/, N_Vector /*tmp2B*/,
                   N_Vector /*tmp3B*/) {

    auto model = static_cast<Model_DAE *>(user_data);
    model->fJB(t, cj, x, dx, xB, dxB, JB);
    return model->checkFinite(gsl::make_span(JB), "Jacobian");
}

/**
 * @brief J in sparse form (for sparse solvers from the SuiteSparse Package)
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
int fJSparse(realtype t, realtype cj, N_Vector x, N_Vector dx,
                        N_Vector /*xdot*/, SUNMatrix J, void *user_data,
                        N_Vector /*tmp1*/, N_Vector /*tmp2*/,
                        N_Vector /*tmp3*/) {
    auto model = static_cast<Model_DAE *>(user_data);
    model->fJSparse(t, cj, x, dx, J);
    return model->checkFinite(gsl::make_span(J), "Jacobian");
}

/**
 * @brief JB in sparse form (for sparse solvers from the SuiteSparse Package)
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
int fJSparseB(realtype t, realtype cj, N_Vector x, N_Vector dx,
                         N_Vector xB, N_Vector dxB, N_Vector /*xBdot*/,
                         SUNMatrix JB, void *user_data, N_Vector /*tmp1B*/,
                         N_Vector /*tmp2B*/, N_Vector /*tmp3B*/) {
    auto model = static_cast<Model_DAE *>(user_data);
    model->fJSparseB(t, cj, x, dx, xB, dxB, JB);
    return model->checkFinite(gsl::make_span(JB), "Jacobian");
}

/**
 * @brief J in banded form (for banded solvers)
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
int fJBand(realtype t, realtype cj, N_Vector x, N_Vector dx,
                      N_Vector xdot, SUNMatrix J, void *user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    return fJ(t, cj, x, dx, xdot, J, user_data, tmp1, tmp2, tmp3);
}

/**
 * @brief JB in banded form (for banded solvers)
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
int fJBandB(realtype t, realtype cj, N_Vector x, N_Vector dx,
                       N_Vector xB, N_Vector dxB, N_Vector xBdot, SUNMatrix JB,
                       void *user_data, N_Vector tmp1B, N_Vector tmp2B,
                       N_Vector tmp3B) {
    return fJB(t, cj, x, dx, xB, dxB, xBdot, JB, user_data, tmp1B, tmp2B,
               tmp3B);
}

/**
 * @brief Matrix vector product of J with a vector v (for iterative solvers)
 * @param t timepoint @type realtype
 * @param cj scaling factor, inverse of the step size
 * @param x Vector with the states
 * @param dx Vector with the derivative states
 * @param xdot Vector with the right hand side
 * @param v Vector with which the Jacobian is multiplied
 * @param Jv Vector to which the Jacobian vector product will be written
 * @param user_data object with user input @type Model_DAE
 * @param tmp1 temporary storage vector
 * @param tmp2 temporary storage vector
 * @return status flag indicating successful execution
 **/
int fJv(realtype t, N_Vector x, N_Vector dx, N_Vector /*xdot*/,
                   N_Vector v, N_Vector Jv, realtype cj, void *user_data,
                   N_Vector /*tmp1*/, N_Vector /*tmp2*/) {

    auto model = static_cast<Model_DAE *>(user_data);
    model->fJv(t, x, dx, v, Jv, cj);
    return model->checkFinite(gsl::make_span(Jv), "Jacobian");
}

/**
 * @brief Matrix vector product of JB with a vector v (for iterative solvers)
 * @param t timepoint
 * @param x Vector with the states
 * @param dx Vector with the derivative states
 * @param xB Vector with the adjoint states
 * @param dxB Vector with the adjoint derivative states
 * @param xBdot Vector with the adjoint right hand side
 * @param vB Vector with which the Jacobian is multiplied
 * @param JvB Vector to which the Jacobian vector product will be written
 * @param cj scalar in Jacobian (inverse stepsize)
 * @param user_data object with user input @type Model_DAE
 * @param tmpB1 temporary storage vector
 * @param tmpB2 temporary storage vector
 * @return status flag indicating successful execution
 **/
int fJvB(realtype t, N_Vector x, N_Vector dx, N_Vector xB,
                    N_Vector dxB, N_Vector /*xBdot*/, N_Vector vB, N_Vector JvB,
                    realtype cj, void *user_data, N_Vector /*tmpB1*/,
                    N_Vector /*tmpB2*/) {

    auto model = static_cast<Model_DAE *>(user_data);
    model->fJvB(t, x, dx, xB, dxB, vB, JvB, cj);
    return model->checkFinite(gsl::make_span(JvB), "Jacobian");
}

/**
 * @brief Event trigger function for events
 * @param t timepoint
 * @param x Vector with the states
 * @param dx Vector with the derivative states
 * @param root array with root function values
 * @param user_data object with user input @type Model_DAE
 * @return status flag indicating successful execution
 */
int froot(realtype t, N_Vector x, N_Vector dx, realtype *root,
                     void *user_data) {
    auto model = static_cast<Model_DAE *>(user_data);
    model->froot(t, x, dx, gsl::make_span<realtype>(root, model->ne));
    return model->checkFinite(gsl::make_span<realtype>(root, model->ne),
                              "root function");
}

/**
 * @brief Residual function of the DAE
 * @param t timepoint
 * @param x Vector with the states
 * @param dx Vector with the derivative states
 * @param xdot Vector with the right hand side
 * @param user_data object with user input @type Model_DAE
 * @return status flag indicating successful execution
 */
int fxdot(realtype t, N_Vector x, N_Vector dx, N_Vector xdot,
                     void *user_data) {
    auto model = static_cast<Model_DAE *>(user_data);

    if (t > 1e200 && !model->app->checkFinite(gsl::make_span(x), "fxdot")) {
        /* when t is large (typically ~1e300), CVODES may pass all NaN x
           to fxdot from which we typically cannot recover. To save time
           on normal execution, we do not always want to check finiteness
           of x, but only do so when t is large and we expect problems. */
        return AMICI_UNRECOVERABLE_ERROR;
    }

    model->fxdot(t, x, dx, xdot);
    return model->checkFinite(gsl::make_span(xdot), "fxdot");
}

/**
 * @brief Right hand side of differential equation for adjoint state xB
 * @param t timepoint
 * @param x Vector with the states
 * @param dx Vector with the derivative states
 * @param xB Vector with the adjoint states
 * @param dxB Vector with the adjoint derivative states
 * @param xBdot Vector with the adjoint right hand side
 * @param user_data object with user input @type Model_DAE
 * @return status flag indicating successful execution
 */
int fxBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB,
           N_Vector dxB, N_Vector xBdot, void *user_data) {

    auto model = static_cast<Model_DAE *>(user_data);
    model->fxBdot(t, x, dx, xB, dxB, xBdot);
    return model->checkFinite(gsl::make_span(xBdot), "xBdot");
}

/**
 * @brief Right hand side of integral equation for quadrature states qB
 * @param t timepoint
 * @param x Vector with the states
 * @param dx Vector with the derivative states
 * @param xB Vector with the adjoint states
 * @param dxB Vector with the adjoint derivative states
 * @param qBdot Vector with the adjoint quadrature right hand side
 * @param user_data pointer to temp data object @type Model_DAE
 * @return status flag indicating successful execution
 */
int fqBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB,
                      N_Vector dxB, N_Vector qBdot, void *user_data) {

    auto model = static_cast<Model_DAE *>(user_data);
    model->fqBdot(t, x, dx, xB, dxB, qBdot);
    return model->checkFinite(gsl::make_span(qBdot), "qBdot");

}


/**
 * @brief Right hand side of differential equation for adjoint state xB
 * when simulating in steadystate mode
 * @param t timepoint
 * @param xB Vector with the adjoint states
 * @param dxB Vector with the adjoint derivative states
 * @param xBdot Vector with the adjoint right hand side
 * @param user_data object with user input @type Model_DAE
 * @return status flag indicating successful execution
 */
static int fxBdot_ss(realtype t, N_Vector xB, N_Vector dxB, N_Vector xBdot,
                     void *user_data) {
    auto model = static_cast<Model_DAE *>(user_data);
    model->fxBdot_ss(t, xB, dxB, xBdot);
    return model->checkFinite(gsl::make_span(xBdot), "xBdot_ss");
}


/**
 * @brief Right hand side of integral equation for quadrature states qB
 * when simulating in steadystate mode
 * @param t timepoint
 * @param xB Vector with the adjoint states
 * @param dxB Vector with the adjoint derivative states
 * @param qBdot Vector with the adjoint quadrature right hand side
 * @param user_data pointer to temp data object
 * @return status flag indicating successful execution
 */
static int fqBdot_ss(realtype t, N_Vector xB, N_Vector dxB, N_Vector qBdot,
                     void *user_data) {
    auto model = static_cast<Model_DAE *>(user_data);
    model->fqBdot_ss(t, xB, dxB, qBdot);
    return model->checkFinite(gsl::make_span(qBdot), "qBdot_ss");
}

/**
 * @brief JB in sparse form for steady state case
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
    static int fJSparseB_ss(realtype /*t*/, realtype /*cj*/, N_Vector /*x*/,
                            N_Vector /*dx*/, N_Vector xBdot, SUNMatrix JB,
                            void *user_data, N_Vector /*tmp1*/,
                            N_Vector /*tmp2*/, N_Vector /*tmp3*/) {
    auto model = static_cast<Model_DAE *>(user_data);
    model->fJSparseB_ss(JB);
    return model->checkFinite(gsl::make_span(xBdot), "JSparseB_ss");
}

/**
 * @brief Right hand side of differential equation for state sensitivities sx
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
int fsxdot(int /*Ns*/, realtype t, N_Vector x, N_Vector dx,
                      N_Vector /*xdot*/, N_Vector *sx, N_Vector *sdx,
                      N_Vector *sxdot, void *user_data, N_Vector /*tmp1*/,
                      N_Vector /*tmp2*/, N_Vector /*tmp3*/) {

    auto model = static_cast<Model_DAE *>(user_data);

    for (int ip = 0; ip < model->nplist(); ip++) {
        model->fsxdot(t, x, dx, ip, sx[ip], sdx[ip], sxdot[ip]);
        if (model->checkFinite(gsl::make_span(sxdot[ip]), "sxdot")
                != AMICI_SUCCESS)
            return AMICI_RECOVERABLE_ERROR;
    }

    return AMICI_SUCCESS;
}

} // namespace amici
