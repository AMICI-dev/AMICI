#ifndef AMICI_SOLVER_IDAS_h
#define AMICI_SOLVER_IDAS_h

#include "amici/solver.h"

#include <sundials/sundials_matrix.h>
#include <idas/idas_impl.h>

namespace amici {

class IDASolver : public Solver {
  public:
    IDASolver() = default;
    ~IDASolver() override = default;

    /**
     * @brief Clone this instance
     * @return The clone
     */
    virtual Solver* clone() const override;

    void reInitPostProcessF(realtype *t, AmiVector *yout, AmiVector *ypout,
                            realtype tnext) override;
    
    void reInitPostProcessB(int which, realtype *t, AmiVector *yBout,
                            AmiVector *ypBout, realtype tnext) override;
    
    void reInit(realtype t0, AmiVector *yy0, AmiVector *yp0) override;

    void sensReInit(AmiVectorArray *yS0, AmiVectorArray *ypS0) override;

    void reInitB(int which, realtype tB0, AmiVector *yyB0,
                   AmiVector *ypB0) override;

    void quadReInitB(int which, AmiVector *yQB0) override;

    void quadSStolerancesB(int which, realtype reltolQB,
                             realtype abstolQB) override;

    int solve(realtype tout, AmiVector *yret, AmiVector *ypret, realtype *tret,
                 int itask) override;

    int solveF(realtype tout, AmiVector *yret, AmiVector *ypret, realtype *tret,
                  int itask, int *ncheckPtr) override;

    void solveB(realtype tBout, int itaskB) override;

    void getRootInfo(int *rootsfound) const override;

    void getDky(realtype t, int k, AmiVector *dky) const override;

    void getSens(realtype *tret, AmiVectorArray *yySout) const override;

    void getB(int which, realtype *tret, AmiVector *yy, AmiVector *yp) const override;

    void getQuadB(int which, realtype *tret, AmiVector *qB) const override;

    void calcIC(realtype tout1, AmiVector *x, AmiVector *dx) override;

    void calcICB(int which, realtype tout1, AmiVector *xB,
                 AmiVector *dxB) override;

    void setStopTime(realtype tstop) override;

    void turnOffRootFinding() override;

    int nplist() const override;

    int nx() const override;

    const Model *getModel() const override;

    bool getMallocDone() const override;

    bool getAdjMallocDone() const override;

    static int fxdot(realtype t, N_Vector x, N_Vector dx, N_Vector xdot,
                     void *user_data);

    static int fJ(realtype t, realtype cj, N_Vector x, N_Vector dx,
                  N_Vector xdot, SUNMatrix J, void *user_data, N_Vector tmp1,
                  N_Vector tmp2, N_Vector tmp3);

    static int fJSparse(realtype t, realtype cj, N_Vector x, N_Vector dx,
                        N_Vector xdot, SUNMatrix J, void *user_data,
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

    void setLinearSolver() override;

    void setLinearSolverB(int which) override;

    void setNonLinearSolver() override;

    void setNonLinearSolverSens() override;

    void setNonLinearSolverB(int which) override;


  protected:
    
    void reInitPostProcess(void *ami_mem, realtype *t, AmiVector *yout,
                           AmiVector *ypout, realtype tout);

    void allocateSolver() override;

    void setSStolerances(double rtol, double atol) override;

    void setSensSStolerances(double rtol, double *atol) override;

    void setSensErrCon(bool error_corr) override;

    void setQuadErrConB(int which, bool flag) override;

    void setErrHandlerFn() override;

    void setUserData(Model *model) override;

    void setUserDataB(int which, Model *model) override;

    void setMaxNumSteps(long int mxsteps) override;

    void setStabLimDet(int stldet) override;

    void setStabLimDetB(int which, int stldet) override;

    void setId(Model *model) override;

    void setSuppressAlg(bool flag) override;
    
     void resetState(void *ida_mem, N_Vector yy0, N_Vector yp0);

    void setSensParams(realtype *p, realtype *pbar, int *plist) override;

    void adjInit() override;

    void allocateSolverB(int *which) override;

    void setMaxNumStepsB(int which, long int mxstepsB) override;

    void setSStolerancesB(int which, realtype relTolB,
                          realtype absTolB) override;

    void diag() override;

    void diagB(int which) override;

    void getNumSteps(void *ami_mem, long int *numsteps) const override;

    void getNumRhsEvals(void *ami_mem, long int *numrhsevals) const override;

    void getNumErrTestFails(void *ami_mem,
                            long int *numerrtestfails) const override;

    void getNumNonlinSolvConvFails(void *ami_mem,
                                   long int *numnonlinsolvconvfails) const override;

    void getLastOrder(void *ami_mem, int *order) const override;

    void *getAdjBmem(void *ami_mem, int which) override;

    void init(AmiVector *x, AmiVector *dx, realtype t) override;

    void binit(int which, AmiVector *xB, AmiVector *dxB, realtype t) override;

    void qbinit(int which, AmiVector *qBdot) override;

    void rootInit(int ne) override;

    void sensInit1(AmiVectorArray *sx, AmiVectorArray *sdx, int nplist) override;

    void setDenseJacFn() override;

    void setSparseJacFn() override;

    void setBandJacFn() override;

    void setJacTimesVecFn() override;

    void setDenseJacFnB(int which) override;

    void setSparseJacFnB(int which) override;

    void setBandJacFnB(int which) override;

    void setJacTimesVecFnB(int which) override;

    static int fJB(realtype t, realtype cj, N_Vector x,
                   N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector xBdot,
                   SUNMatrix JB, void *user_data, N_Vector tmp1B,
                   N_Vector tmp2B, N_Vector tmp3B);

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

    static int fqBdot(realtype t, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector qBdot,
                      void *user_data);

    static int fsxdot(int Ns, realtype t, N_Vector x, N_Vector dx, N_Vector xdot,
                      N_Vector *sx, N_Vector *sdx, N_Vector *sxdot, void *user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
};

} // namespace amici

#endif /* idawrap_h */
