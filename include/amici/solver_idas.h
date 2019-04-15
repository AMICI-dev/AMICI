#ifndef AMICI_SOLVER_IDAS_h
#define AMICI_SOLVER_IDAS_h

#include "amici/solver.h"

#include <sundials/sundials_matrix.h>

namespace amici {
class ExpData;
class ReturnData;
class Model_DAE;
class IDASolver;
} // namespace amici

// for serialization friend in Solver
namespace boost {
namespace serialization {
template <class Archive>
void serialize(Archive &ar, amici::IDASolver &u, const unsigned int version);
}
} // namespace boost::serialization

namespace amici {

class IDASolver : public Solver {
  public:
    IDASolver() = default;
    ~IDASolver() override = default;

    /**
     * @brief Clone this instance
     * @return The clone
     */
    virtual Solver *clone() const override;

    void reInitPostProcessF(const realtype tnext) const override;

    void reInitPostProcessB(const realtype tnext) const override;

    void reInit(const realtype t0, const AmiVector &yy0,
                const AmiVector &yp0) const override;
    
    void sensReInit(const AmiVectorArray &yyS0,
                    const AmiVectorArray &ypS0) const override;
    
    void reInitB(const int which, const realtype tB0,
                 const AmiVector &yyB0, const AmiVector &ypB0) const override;
    
    void quadReInitB(const int which, const AmiVector &yQB0) const override;

    void quadSStolerancesB(const int which, const realtype reltolQB,
                           const realtype abstolQB) const override;

    int solve(const realtype tout, const int itask) const override;

    int solveF(const realtype tout, const int itask,
               int *ncheckPtr) const override;

    void solveB(const realtype tBout, const int itaskB) const override;

    void getRootInfo(int *rootsfound) const override;

    void getDky(const realtype t, int k) const override;

    void getSens() const override;

    void getSensDky(const realtype t, int k) const override;

    void getB(int which) const override;

    void getDkyB(const realtype t, int k, int which) const override;

    void getQuadB(int which) const override;

    void getQuadDkyB(const realtype t, int k, int which) const override;

    void calcIC(realtype tout1) const override;

    void calcICB(int which, realtype tout1) const override;

    void setStopTime(realtype tstop) const override;

    void turnOffRootFinding() const override;

    const Model *getModel() const override;

    static int fxdot(realtype t, N_Vector x, N_Vector dx, N_Vector xdot,
                     void *user_data);

    static int fJ(realtype t, realtype cj, N_Vector x, N_Vector dx,
                  N_Vector xdot, SUNMatrix J, void *user_data, N_Vector tmp1,
                  N_Vector tmp2, N_Vector tmp3);

    static int fJSparse(realtype t, realtype cj, N_Vector x, N_Vector dx,
                        N_Vector xdot, SUNMatrix J, void *user_data,
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

    void setLinearSolver() const override;

    void setLinearSolverB(int which) const override;

    void setNonLinearSolver() const override;

    void setNonLinearSolverSens() const override;

    void setNonLinearSolverB(int which) const override;

  protected:
    void reInitPostProcess(void *ami_mem, realtype *t, AmiVector *yout,
                           AmiVector *ypout, realtype tout) const;

    void allocateSolver() const override;

    void setSStolerances(const realtype rtol,
                         const realtype atol) const override;

    void setSensSStolerances(const realtype rtol,
                             const realtype *atol) const override;

    void setSensErrCon(const bool error_corr) const override;

    void setQuadErrConB(const int which, const bool flag) const override;

    void setErrHandlerFn() const override;

    void setUserData(Model *model) const override;

    void setUserDataB(const int which, Model *model) const override;

    void setMaxNumSteps(const long int mxsteps) const override;

    void setStabLimDet(const int stldet) const override;

    void setStabLimDetB(const int which, const int stldet) const override;

    void setId(const Model *model) const override;

    void setSuppressAlg(bool flag) const override;

    void resetState(void *ida_mem, const N_Vector yy0,
                    const N_Vector yp0) const;

    void setSensParams(const realtype *p, const realtype *pbar,
                       const int *plist) const override;

    void adjInit() const override;

    void allocateSolverB(int *which) const override;

    void setMaxNumStepsB(const int which,
                         const long int mxstepsB) const override;

    void setSStolerancesB(const int which, const realtype relTolB,
                          const realtype absTolB) const override;

    void diag() const override;

    void diagB(const int which) const override;

    void getNumSteps(const void *ami_mem, long int *numsteps) const override;

    void getNumRhsEvals(const void *ami_mem,
                        long int *numrhsevals) const override;

    void getNumErrTestFails(const void *ami_mem,
                            long int *numerrtestfails) const override;

    void
    getNumNonlinSolvConvFails(const void *ami_mem,
                              long int *numnonlinsolvconvfails) const override;

    void getLastOrder(const void *ami_mem, int *order) const override;

    void *getAdjBmem(void *ami_mem, int which) const override;

    void init(const realtype t0, const AmiVector &x0,
              const AmiVector &dx0) const override;

    void sensInit1(const AmiVectorArray &sx0, const AmiVectorArray &sdx0) const override;

    void binit(const int which, const realtype tf,
               const AmiVector &xB0, const AmiVector &dxB0) const override;

    void qbinit(const int which, const AmiVector &xQB0) const override;

    void rootInit(int ne) const override;

    void setDenseJacFn() const override;

    void setSparseJacFn() const override;

    void setBandJacFn() const override;

    void setJacTimesVecFn() const override;

    void setDenseJacFnB(const int which) const override;

    void setSparseJacFnB(const int which) const override;

    void setBandJacFnB(const int which) const override;

    void setJacTimesVecFnB(const int which) const override;

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

    static int fsxdot(int Ns, realtype t, N_Vector x, N_Vector dx,
                      N_Vector xdot, N_Vector *sx, N_Vector *sdx,
                      N_Vector *sxdot, void *user_data, N_Vector tmp1,
                      N_Vector tmp2, N_Vector tmp3);
};

} // namespace amici

#endif /* idawrap_h */
