#ifndef AMICI_SOLVER_CVODES_h
#define AMICI_SOLVER_CVODES_h

#include "include/amici_solver.h"
#include "include/amici_exception.h"
#include "include/amici_defines.h"
#include <cvodes/cvodes_dense.h>
#include <sundials/sundials_sparse.h>

namespace amici {

class UserData;
class ExpData;
class ReturnData;
class Model_ODE;

class CVodeSolver : public Solver {
  public:
    CVodeSolver(const Model_ODE *model){
        this->model = model;
    }

    void *AMICreate(int lmm, int iter) override;

    void AMISStolerances(double rtol, double atol) override;

    void AMISensEEtolerances() override;

    void AMISetSensErrCon(bool error_corr) override;

    void AMISetQuadErrConB(int which, bool flag) override;

    void AMIGetRootInfo(int *rootsfound) override;

    void AMISetErrHandlerFn() override;

    void AMISetUserData(void *user_data) override;

    void AMISetUserDataB(int which, void *user_data) override;

    void AMISetMaxNumSteps(long int mxsteps) override;

    void AMISetStabLimDet(int stldet) override;

    void AMISetStabLimDetB(int which, int stldet) override;

    void AMISetId(Model *model) override;

    void AMISetSuppressAlg(bool flag) override;

    void AMIReInit(realtype t0, AmiVector yy0, AmiVector yp0) override;

    void AMISensReInit(int ism, AmiVectorArray yS0, AmiVectorArray ypS0) override;

    void AMISetSensParams(realtype *p, realtype *pbar, int *plist) override;

    void AMIGetDky(realtype t, int k, AmiVector dky) override;

    void AMIGetSens(realtype *tret, AmiVectorArray yySout) override;

    void AMIFree() override;

    void AMIAdjInit(long int steps, int interp) override;

    void AMICreateB(int lmm, int iter, int *which) override;

    void AMIReInitB(int which, realtype tB0, AmiVector yyB0,
                   AmiVector ypB0) override;

    void AMISStolerancesB(int which, realtype relTolB,
                         realtype absTolB) override;

    void AMIQuadReInitB(int which, AmiVector yQB0) override;

    void AMIQuadSStolerancesB(int which, realtype reltolQB,
                             realtype abstolQB) override;

    int AMISolve(realtype tout, AmiVector yret, AmiVector ypret, realtype *tret,
                 int itask) override;

    int AMISolveF(realtype tout, AmiVector yret, AmiVector ypret, realtype *tret,
                  int itask, int *ncheckPtr) override;

    void AMISolveB(realtype tBout, int itaskB) override;

    void AMISetMaxNumStepsB(int which, long int mxstepsB) override;

    void AMIGetB(int which, realtype *tret, AmiVector yy, AmiVector yp) override;

    void AMIGetQuadB(int which, realtype *tret, AmiVector qB) override;

    void AMIDense(int nx) override;

    void AMIDenseB(int which, int nx) override;

    void AMIBand(int nx, int ubw, int lbw) override;

    void AMIBandB(int which, int nx, int ubw, int lbw) override;

    void AMIDiag() override;

    void AMIDiagB(int which) override;

    void AMISpgmr(int prectype, int maxl) override;

    void AMISpgmrB(int which, int prectype, int maxl) override;

    void AMISpbcg(int prectype, int maxl) override;

    void AMISpbcgB(int which, int prectype, int maxl) override;

    void AMISptfqmr(int prectype, int maxl) override;

    void AMISptfqmrB(int which, int prectype, int maxl) override;

    void AMIKLU(int nx, int nnz, int sparsetype) override;

    void AMIKLUSetOrdering(int ordering) override;

    void AMIKLUSetOrderingB(int which, int ordering) override;

    void AMIKLUB(int which, int nx, int nnz, int sparsetype) override;

    void AMIGetNumSteps(void *ami_mem, long int *numsteps) override;

    void AMIGetNumRhsEvals(void *ami_mem, long int *numrhsevals) override;

    void AMIGetNumErrTestFails(void *ami_mem,
                              long int *numerrtestfails) override;

    void AMIGetNumNonlinSolvConvFails(void *ami_mem,
                                     long int *numnonlinsolvconvfails) override;

    void AMIGetLastOrder(void *ami_ami_mem, int *order) override;

    void *AMIGetAdjBmem(void *ami_mem, int which) override;


    void AMICalcIC(realtype tout1, AmiVector x, AmiVector dx) override;

    void AMICalcICB(int which, realtype tout1, AmiVector xB,
                   AmiVector dxB) override;

    void AMISetStopTime(realtype tstop) override;

    void turnOffRootFinding() override;

    ~CVodeSolver();

  protected:
    
    const Model_ODE *model;
    
    void init(AmiVector x, AmiVector dx, realtype t) override;

    void binit(int which, AmiVector xB, AmiVector dxB, realtype t) override;

    void qbinit(int which, AmiVector qBdot) override;

    void rootInit(int ne) override;

    void sensInit1(AmiVectorArray sx, AmiVectorArray sdx, const UserData *udata) override;

    void setDenseJacFn() override;

    void setSparseJacFn() override;

    void setBandJacFn() override;

    void setJacTimesVecFn() override;

    void setDenseJacFnB(int which) override;

    void setSparseJacFnB(int which) override;

    void setBandJacFnB(int which) override;

    void setJacTimesVecFnB(int which) override;
};

} // namespace amici

#endif /* CVodewrap_h */
