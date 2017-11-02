#ifndef idawrap_h
#define idawrap_h

#include "include/amici_solver.h"
#include <cvodes/cvodes_dense.h>
#include <sundials/sundials_sparse.h>

namespace amici {

class IDASolver : public Solver {
  public:
    IDASolver();

    void *AMICreate(int lmm, int iter) override;

    int AMISStolerances(double rtol, double atol) override;

    int AMISensEEtolerances() override;

    int AMISetSensErrCon(bool error_corr) override;

    int AMISetQuadErrConB(int which, bool flag) override;

    int AMIGetRootInfo(int *rootsfound) override;

    int AMISetErrHandlerFn() override;

    int AMISetUserData(void *user_data) override;

    int AMISetUserDataB(int which, void *user_data) override;

    int AMISetMaxNumSteps(long int mxsteps) override;

    int AMISetStabLimDet(int stldet) override;

    int AMISetStabLimDetB(int which, int stldet) override;

    int AMISetId(Model *model) override;

    int AMISetSuppressAlg(bool flag) override;

    int AMIReInit(realtype t0, N_Vector yy0, N_Vector yp0) override;

    int AMISensReInit(int ism, N_Vector *yS0, N_Vector *ypS0) override;

    int AMISetSensParams(realtype *p, realtype *pbar, int *plist) override;

    int AMIGetDky(realtype t, int k, N_Vector dky) override;

    int AMIGetSens(realtype *tret, N_Vector *yySout) override;

    void AMIFree() override;

    int AMIAdjInit(long int steps, int interp) override;

    int AMICreateB(int lmm, int iter, int *which) override;

    int AMIReInitB(int which, realtype tB0, N_Vector yyB0,
                   N_Vector ypB0) override;

    int AMISStolerancesB(int which, realtype relTolB,
                         realtype absTolB) override;

    int AMIQuadReInitB(int which, N_Vector yQB0) override;

    int AMIQuadSStolerancesB(int which, realtype reltolQB,
                             realtype abstolQB) override;

    int AMISolve(realtype tout, N_Vector yret, N_Vector ypret, realtype *tret,
                 int itask) override;

    int AMISolveF(realtype tout, N_Vector yret, N_Vector ypret, realtype *tret,
                  int itask, int *ncheckPtr) override;

    int AMISolveB(realtype tBout, int itaskB) override;

    int AMISetMaxNumStepsB(int which, long int mxstepsB) override;

    int AMIGetB(int which, realtype *tret, N_Vector yy, N_Vector yp) override;

    int AMIGetQuadB(int which, realtype *tret, N_Vector qB) override;

    int AMIDense(int nx) override;

    int AMIDenseB(int which, int nx) override;

    int AMIBand(int nx, int ubw, int lbw) override;

    int AMIBandB(int which, int nx, int ubw, int lbw) override;

    int AMIDiag() override;

    int AMIDiagB(int which) override;

    int AMISpgmr(int prectype, int maxl) override;

    int AMISpgmrB(int which, int prectype, int maxl) override;

    int AMISpbcg(int prectype, int maxl) override;

    int AMISpbcgB(int which, int prectype, int maxl) override;

    int AMISptfqmr(int prectype, int maxl) override;

    int AMISptfqmrB(int which, int prectype, int maxl) override;

    int AMIKLU(int nx, int nnz, int sparsetype) override;

    int AMIKLUSetOrdering(int ordering) override;

    int AMIKLUSetOrderingB(int which, int ordering) override;

    int AMIKLUB(int which, int nx, int nnz, int sparsetype) override;

    int AMIGetNumSteps(void *ami_mem, long int *numsteps) override;

    int AMIGetNumRhsEvals(void *ami_mem, long int *numrhsevals) override;

    int AMIGetNumErrTestFails(void *ami_mem,
                              long int *numerrtestfails) override;

    int AMIGetNumNonlinSolvConvFails(void *ami_mem,
                                     long int *numnonlinsolvconvfails) override;

    int AMIGetLastOrder(void *ami_mem, int *order) override;

    void *AMIGetAdjBmem(void *ami_mem, int which) override;

    int AMICalcIC(realtype tout1, TempData *tdata) override;

    int AMICalcICB(int which, realtype tout1, N_Vector xB,
                   N_Vector dxB) override;

    int AMISetStopTime(realtype tstop) override;

    int turnOffRootFinding() override;

    static int resultFunction(realtype tt, N_Vector yy, N_Vector yp,
                              N_Vector rr, void *user_data);

    static int resultFunctionB(realtype tt, N_Vector yy, N_Vector yp,
                               N_Vector yyB, N_Vector ypB, N_Vector rrB,
                               void *user_dataB);

    static int rootFunction(realtype t, N_Vector y, N_Vector yp, realtype *gout,
                            void *user_data);

    static int J(long int N, realtype t, realtype c_j, N_Vector y, N_Vector yp,
                 N_Vector r, DlsMat Jac, void *user_data, N_Vector tmp1,
                 N_Vector tmp2, N_Vector tmp3);

    static int fqBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB,
                      N_Vector dxB, N_Vector qBdot, void *user_data);

    static int fsxdot(int Ns, realtype t, N_Vector x, N_Vector xdot,
                      N_Vector dx, N_Vector *sx, N_Vector *sxdot, N_Vector *sdx,
                      void *user_data, N_Vector tmp1, N_Vector tmp2,
                      N_Vector tmp3);

    static int fJSparse(realtype t, realtype cj, N_Vector x, N_Vector dx,
                        N_Vector xdot, SlsMat J, void *user_data, N_Vector tmp1,
                        N_Vector tmp2, N_Vector tmp3);

    static int fJBand(long int N, long int mupper, long int mlower, realtype t,
                      realtype cj, N_Vector x, N_Vector dx, N_Vector xdot,
                      DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2,
                      N_Vector tmp3);

    static int fJv(realtype t, N_Vector x, N_Vector dx, N_Vector xdot,
                   N_Vector v, N_Vector Jv, realtype cj, void *user_data,
                   N_Vector tmp1, N_Vector tmp2);

    static int fJB(long int NeqBdot, realtype t, realtype cj, N_Vector x,
                   N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot,
                   DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B,
                   N_Vector tmp3B);

    static int fJSparseB(realtype t, realtype cj, N_Vector x, N_Vector dx,
                         N_Vector xB, N_Vector dxB, N_Vector xBdot, SlsMat JB,
                         void *user_data, N_Vector tmp1B, N_Vector tmp2B,
                         N_Vector tmp3B);

    static int fJBandB(long int NeqBdot, long int mupper, long int mlower,
                       realtype t, realtype cj, N_Vector x, N_Vector dx,
                       N_Vector xB, N_Vector dxB, N_Vector xBdot, DlsMat JB,
                       void *user_data, N_Vector tmp1B, N_Vector tmp2B,
                       N_Vector tmp3B);

    static int fJvB(realtype t, N_Vector x, N_Vector dx, N_Vector xB,
                    N_Vector dxB, N_Vector xBdot, N_Vector vB, N_Vector JvB,
                    realtype cj, void *user_data, N_Vector tmpB1,
                    N_Vector tmpB2);

    ~IDASolver();

  protected:
    int init(N_Vector x, N_Vector dx, realtype t) override;

    int binit(int which, N_Vector xB, N_Vector dxB, realtype t) override;

    int qbinit(int which, N_Vector qBdot) override;

    int rootInit(int ne) override;

    int sensInit1(N_Vector *sx, N_Vector *sdx, const UserData *udata) override;

    int setDenseJacFn() override;

    int setSparseJacFn() override;

    int setBandJacFn() override;

    int setJacTimesVecFn() override;

    int setDenseJacFnB(int which) override;

    int setSparseJacFnB(int which) override;

    int setBandJacFnB(int which) override;

    int setJacTimesVecFnB(int which) override;
};

} // namespace amici

#endif /* idawrap_h */
