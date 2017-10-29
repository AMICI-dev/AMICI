#ifndef AMICI_SOLVER_CVODES_h
#define AMICI_SOLVER_CVODES_h

#include "include/amici_solver.h"
#include "include/amici_exception.h"
#include <cvodes/cvodes_dense.h>
#include <sundials/sundials_sparse.h>

namespace amici {

class UserData;

class CVodeSolver : public Solver {
  public:
    CVodeSolver() = default;

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

    void AMIReInit(realtype t0, N_Vector yy0, N_Vector yp0) override;

    void AMISensReInit(int ism, N_Vector *yS0, N_Vector *ypS0) override;

    void AMISetSensParams(realtype *p, realtype *pbar, int *plist) override;

    void AMIGetDky(realtype t, int k, N_Vector dky) override;

    void AMIGetSens(realtype *tret, N_Vector *yySout) override;

    void AMIFree() override;

    void AMIAdjInit(long int steps, int interp) override;

    void AMICreateB(int lmm, int iter, int *which) override;

    void AMIReInitB(int which, realtype tB0, N_Vector yyB0,
                   N_Vector ypB0) override;

    void AMISStolerancesB(int which, realtype relTolB,
                         realtype absTolB) override;

    void AMIQuadReInitB(int which, N_Vector yQB0) override;

    void AMIQuadSStolerancesB(int which, realtype reltolQB,
                             realtype abstolQB) override;

    int AMISolve(realtype tout, N_Vector yret, N_Vector ypret, realtype *tret,
                 int itask) override;

    int AMISolveF(realtype tout, N_Vector yret, N_Vector ypret, realtype *tret,
                  int itask, int *ncheckPtr) override;

    void AMISolveB(realtype tBout, int itaskB) override;

    void AMISetMaxNumStepsB(int which, long int mxstepsB) override;

    void AMIGetB(int which, realtype *tret, N_Vector yy, N_Vector yp) override;

    void AMIGetQuadB(int which, realtype *tret, N_Vector qB) override;

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

    void AMICalcIC(realtype tout1) override;

    void AMICalcICB(int which, realtype tout1, N_Vector xB,
                   N_Vector dxB) override;

    void AMISetStopTime(realtype tstop) override;

    void turnOffRootFinding() override;

    // Static wrapper functions because cannot pass member functions to solver
    // (CVODES-specific signatures)
    static int residualFunction(realtype t, N_Vector y, N_Vector ydot,
                              void *user_data);

    static int residualFunctionB(realtype t, N_Vector y, N_Vector yB,
                               N_Vector yBdot, void *user_dataB);

    static int rootFunction(realtype t, N_Vector x, realtype *root,
                            void *user_data);

    static int J(long int N, realtype t, N_Vector x, N_Vector xdot, DlsMat J,
                 void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

    static int fqBdot(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot,
                      void *user_data);

    static int fsxdot(int Ns, realtype t, N_Vector x, N_Vector xdot, int ip,
                      N_Vector sx, N_Vector sxdot, void *user_data,
                      N_Vector tmp1, N_Vector tmp2);

    static int fJSparse(realtype t, N_Vector x, N_Vector xdot, SlsMat J,
                        void *user_data, N_Vector tmp1, N_Vector tmp2,
                        N_Vector tmp3);

    static int fJBand(long int N, long int mupper, long int mlower, realtype t,
                      N_Vector x, N_Vector xdot, DlsMat J, void *user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

    static int fJv(N_Vector v, N_Vector Jv, realtype t, N_Vector x,
                   N_Vector xdot, void *user_data, N_Vector tmp);

    static int fJB(long NeqBdot, realtype t, N_Vector x, N_Vector xB,
                   N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B,
                   N_Vector tmp2B, N_Vector tmp3B);

    static int fJSparseB(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot,
                         SlsMat JB, void *user_data, N_Vector tmp1B,
                         N_Vector tmp2B, N_Vector tmp3B);

    static int fJBandB(long NeqBdot, long mupper, long mlower, realtype t,
                       N_Vector x, N_Vector xB, N_Vector xBdot, DlsMat JB,
                       void *user_data, N_Vector tmp1B, N_Vector tmp2B,
                       N_Vector tmp3B);

    static int fJvB(N_Vector vB, N_Vector JvB, realtype t, N_Vector x,
                    N_Vector xB, N_Vector xBdot, void *user_data,
                    N_Vector tmpB);

    ~CVodeSolver();

  protected:
    void init(N_Vector x, N_Vector dx, realtype t) override;

    void binit(int which, N_Vector xB, N_Vector dxB, realtype t) override;

    void qbinit(int which, N_Vector qBdot) override;

    void rootInit(int ne) override;

    void sensInit1(N_Vector *sx, N_Vector *sdx, const UserData *udata) override;

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
