#ifndef AMICI_SOLVER_CVODES_h
#define AMICI_SOLVER_CVODES_h

#include "amici/defines.h"
#include "amici/solver.h"
#include "amici/vector.h"

#include <sundials/sundials_matrix.h>

namespace amici {
class ExpData;
class ReturnData;
class Model_ODE;
class CVodeSolver;
} // namespace amici

// for serialization friend in Solver
namespace boost {
namespace serialization {
template <class Archive>
void serialize(Archive &ar, amici::CVodeSolver &u, const unsigned int version);
}
} // namespace boost

namespace amici {

class CVodeSolver : public Solver {
  public:
    CVodeSolver() = default;

    ~CVodeSolver() override = default;

    /**
     * @brief Clone this instance
     * @return The clone
     */
    virtual Solver *clone() const override;
    
    void reInit(const realtype t0, const AmiVector &yy0,
                const AmiVector &yp0) const override;
    
    void sensReInit(const AmiVectorArray &yyS0,
                    const AmiVectorArray &ypS0) const override;
    
    void reInitB(const int which, const realtype tB0,
                 const AmiVector &yyB0, const AmiVector &ypB0) const override;
    
    void quadReInitB(const int which, const AmiVector &yQB0) const override;

    int solve(const realtype tout, const int itask) const override;

    int solveF(const realtype tout, const int itask,
               int *ncheckPtr) const override;

    void solveB(const realtype tBout, int itaskB) const override;

    void getDky(const realtype t, const int k) const override;

    void getSensDky(const realtype t, const int k) const override;

    void getQuadDkyB(const realtype t, const int k,
                     const int which) const override;

    void getDkyB(const realtype t, const int k, const int which) const override;

    void getRootInfo(int *rootsfound) const override;

    void setStopTime(const realtype tstop) const override;

    void turnOffRootFinding() const override;

    const Model *getModel() const override;

    static int fxdot(realtype t, N_Vector x, N_Vector xdot, void *user_data);

    static int fJSparse(realtype t, N_Vector x, N_Vector xdot, SUNMatrix J,
                        void *user_data, N_Vector tmp1, N_Vector tmp2,
                        N_Vector tmp3);

    static int fJ(realtype t, N_Vector x, N_Vector xdot, SUNMatrix J,
                  void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

    using Solver::setLinearSolver;
    using Solver::setLinearSolverB;

    void setLinearSolver() const override;

    void setLinearSolverB(const int which) const override;

    void setNonLinearSolver() const override;

    void setNonLinearSolverSens() const override;

    void setNonLinearSolverB(const int which) const override;

  protected:
    
    void calcIC(const realtype tout1) const override;

    void calcICB(const int which, const realtype tout1) const override;

    void getB(const int which) const override;

    void getSens() const override;

    void getQuadB(const int which) const override;

    void reInitPostProcessF(const realtype tnext) const override;

    void reInitPostProcessB(const realtype tnext) const override;

    void reInitPostProcess(void *ami_mem, realtype *t, AmiVector *yout,
                           realtype tout) const;

    void allocateSolver() const override;

    void setSStolerances(const double rtol, const double atol) const override;

    void setSensSStolerances(const double rtol,
                             const double *atol) const override;

    void setSensErrCon(const bool error_corr) const override;

    void setQuadErrConB(const int which, const bool flag) const override;

    void setErrHandlerFn() const override;

    void setUserData(Model *model) const override;

    void setUserDataB(const int which, Model *model) const override;

    void setMaxNumSteps(const long int mxsteps) const override;

    void setStabLimDet(const int stldet) const override;

    void setStabLimDetB(const int which, int stldet) const override;

    void setId(const Model *model) const override;

    void setSuppressAlg(const bool flag) const override;

    void resetState(void *cv_mem, const N_Vector y0) const;

    void setSensParams(const realtype *p, const realtype *pbar,
                       const int *plist) const override;

    void adjInit() const override;

    void allocateSolverB(int *which) const override;

    void setSStolerancesB(const int which, const realtype relTolB,
                          const realtype absTolB) const override;

    void quadSStolerancesB(const int which, const realtype reltolQB,
                           realtype abstolQB) const override;

    void setMaxNumStepsB(int which, long int mxstepsB) const override;

    void diag() const override;

    void diagB(int which) const override;

    void getNumSteps(const void *ami_mem, long int *numsteps) const override;

    void getNumRhsEvals(const void *ami_mem,
                        long int *numrhsevals) const override;

    void getNumErrTestFails(const void *ami_mem,
                            long int *numerrtestfails) const override;

    void
    getNumNonlinSolvConvFails(const void *ami_mem,
                              long int *numnonlinsolvconvfails) const override;

    void getLastOrder(const void *ami_ami_mem, int *order) const override;

    void *getAdjBmem(void *ami_mem, int which) const override;

    template <class Archive>
    friend void boost::serialization::serialize(Archive &ar, CVodeSolver &r,
                                                const unsigned int version);

    friend bool operator==(const CVodeSolver &a, const CVodeSolver &b);

    void init(const realtype t0, const AmiVector &x0, const AmiVector &dx0)
    const override;

    void sensInit1(const AmiVectorArray &sx0, const AmiVectorArray &sdx0)
    const override;

    void binit(const int which, const realtype tf, const AmiVector &xB0,
               const AmiVector &dxB0) const override;

    void qbinit(const int which, const AmiVector &xQB0) const override;

    void rootInit(const int ne) const override;

    void setDenseJacFn() const override;

    void setSparseJacFn() const override;

    void setBandJacFn() const override;

    void setJacTimesVecFn() const override;

    void setDenseJacFnB(const int which) const override;

    void setSparseJacFnB(const int which) const override;

    void setBandJacFnB(const int which) const override;

    void setJacTimesVecFnB(const int which) const override;

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

    static int fJDiag(realtype t, N_Vector JDiag, N_Vector x, void *user_data);

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
};

} // namespace amici

#endif /* CVodewrap_h */
