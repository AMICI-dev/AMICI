#ifndef idawrap_h
#define idawrap_h

#include <sundials/sundials_sparse.h>
#include <cvodes/cvodes_dense.h>
#include "include/amici_solver.h"

class IDASolver : public Solver
{
public:
    IDASolver();

    int wrap_init(N_Vector x, N_Vector dx, realtype t);

    int wrap_binit(int which, N_Vector xB, N_Vector dxB, realtype t);

    int wrap_qbinit(int which, N_Vector qBdot);

    int wrap_RootInit(int ne);

    int wrap_SensInit1(N_Vector *sx, N_Vector *sdx, const UserData *udata);

    int wrap_SetDenseJacFn();

    int wrap_SetSparseJacFn();

    int wrap_SetBandJacFn();

    int wrap_SetJacTimesVecFn();

    int wrap_SetDenseJacFnB(int which);

    int wrap_SetSparseJacFnB(int which);

    int wrap_SetBandJacFnB(int which);

    int wrap_SetJacTimesVecFnB(int which);

    void *AMICreate(int lmm, int iter);

    int AMISStolerances(double rtol,double atol);

    int AMISensEEtolerances();

    int AMISetSensErrCon(bool error_corr);

    int AMISetQuadErrConB(int which, bool flag);

    int AMIGetRootInfo(int *rootsfound);

    int AMISetErrHandlerFn();

    int AMISetUserData(void *user_data);

    int AMISetUserDataB(int which, void *user_data);

    int AMISetMaxNumSteps(long int mxsteps);

    int AMISetStabLimDet(int stldet);

    int AMISetStabLimDetB(int which, int stldet);

    int AMISetId(Model *model);

    int AMISetSuppressAlg(bool flag);

    int AMIReInit(realtype t0, N_Vector yy0, N_Vector yp0);

    int AMISensReInit(int ism, N_Vector *yS0, N_Vector *ypS0);

    int AMISetSensParams(realtype *p, realtype *pbar, int *plist);

    int AMIGetDky(realtype t, int k, N_Vector dky);

    int AMIGetSens(realtype *tret, N_Vector *yySout);

    void AMIFree();

    int AMIAdjInit(long int steps, int interp);

    int AMICreateB(int lmm, int iter, int *which);

    int AMIReInitB(int which, realtype tB0, N_Vector yyB0, N_Vector ypB0);

    int AMISStolerancesB(int which, realtype relTolB, realtype absTolB);

    int AMIQuadReInitB(int which, N_Vector yQB0);

    int AMIQuadSStolerancesB(int which, realtype reltolQB, realtype abstolQB);

    int AMISolve(realtype tout, N_Vector yret, N_Vector ypret, realtype *tret, int itask);

    int AMISolveF(realtype tout, N_Vector yret, N_Vector ypret, realtype *tret, int itask, int *ncheckPtr);

    int AMISolveB(realtype tBout, int itaskB);

    int AMISetMaxNumStepsB(int which, long int mxstepsB);

    int AMIGetB(int which, realtype *tret, N_Vector yy, N_Vector yp);

    int AMIGetQuadB(int which, realtype *tret, N_Vector qB);

    int AMIDense(int nx);

    int AMIDenseB(int which, int nx);

    int AMIBand(int nx, int ubw, int lbw);

    int AMIBandB(int which, int nx, int ubw, int lbw);

    int AMIDiag();

    int AMIDiagB(int which);

    int AMISpgmr(int prectype, int maxl);

    int AMISpgmrB(int which, int prectype, int maxl);

    int AMISpbcg(int prectype, int maxl);

    int AMISpbcgB(int which, int prectype, int maxl);

    int AMISptfqmr(int prectype, int maxl);

    int AMISptfqmrB(int which, int prectype, int maxl);

    int AMIKLU(int nx, int nnz, int sparsetype);

    int AMIKLUSetOrdering(int ordering);

    int AMIKLUSetOrderingB(int which, int ordering);

    int AMIKLUB(int which, int nx, int nnz, int sparsetype);

    int AMIGetNumSteps(void *ami_mem, long int *numsteps);

    int AMIGetNumRhsEvals(void *ami_mem, long int *numrhsevals);

    int AMIGetNumErrTestFails(void *ami_mem, long int *numerrtestfails);

    int AMIGetNumNonlinSolvConvFails(void *ami_mem, long int *numnonlinsolvconvfails);

    int AMIGetLastOrder(void *ami_mem, int *order);

    void *AMIGetAdjBmem(void *ami_mem, int which);

    int AMICalcIC(realtype tout1);

    int AMICalcICB(int which, realtype tout1, N_Vector xB, N_Vector dxB);

    int AMISetStopTime(realtype tstop);

    int turnOffRootFinding();


    static int resultFunction(realtype tt, N_Vector yy, N_Vector yp,
                              N_Vector rr, void *user_data);

    static int resultFunctionB(realtype tt,
                               N_Vector yy, N_Vector yp,
                               N_Vector yyB, N_Vector ypB,
                               N_Vector rrB, void *user_dataB);


    static int rootFunction(realtype t, N_Vector y, N_Vector yp,
                             realtype *gout, void *user_data);

    static int J(long int N, realtype t, realtype c_j,
                                    N_Vector y, N_Vector yp, N_Vector r,
                                    DlsMat Jac, void *user_data,
                                    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);


    static int fqBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector qBdot, void *user_data);

    static int fsxdot(int Ns, realtype t,
                      N_Vector x, N_Vector xdot, N_Vector dx,
                      N_Vector *sx, N_Vector *sxdot, N_Vector *sdx, void *user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

    static int fJSparse(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, SlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

    static int fJBand(long int N, long int mupper, long int mlower, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

    static int fJv(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, N_Vector v, N_Vector Jv, realtype cj, void *user_data, N_Vector tmp1, N_Vector tmp2);

    static int fJB(long int NeqBdot, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

    static int fJSparseB(realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, SlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

    static int fJBandB(long int NeqBdot, long int mupper, long int mlower, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

    static int fJvB(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot, N_Vector vB, N_Vector JvB, realtype cj, void *user_data, N_Vector tmpB1, N_Vector tmpB2);

    ~IDASolver();
};

#endif /* idawrap_h */
