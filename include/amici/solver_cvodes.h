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
void serialize(Archive &ar, amici::CVodeSolver &u, unsigned int version);
}
} // namespace boost::serialization

namespace amici {

class CVodeSolver : public Solver {
  public:
    ~CVodeSolver() override = default;

    /**
     * @brief Clone this instance
     * @return The clone
     */
    Solver *clone() const override;

    void reInit(realtype t0, const AmiVector &yy0,
                const AmiVector &yp0) const override;

    void sensReInit(const AmiVectorArray &yyS0,
                    const AmiVectorArray &ypS0) const override;

    void sensToggleOff() const override;

    void reInitB(int which, realtype tB0,
                 const AmiVector &yyB0, const AmiVector &ypB0) const override;

    void quadReInitB(int which, const AmiVector &yQB0) const override;

    int solve(realtype tout, int itask) const override;

    int solveF(realtype tout, int itask,
               int *ncheckPtr_) const override;

    void solveB(realtype tBout, int itaskB) const override;

    void getDky(realtype t_, int k) const override;

    void getSensDky(realtype t_, int k) const override;

    void getQuadDkyB(realtype t_, int k,
                     int which) const override;

    void getDkyB(realtype t_, int k, int which) const override;

    void getRootInfo(int *rootsfound) const override;

    void setStopTime(realtype tstop) const override;

    void turnOffRootFinding() const override;

    const Model *getModel() const override;

    using Solver::setLinearSolver;

    using Solver::setLinearSolverB;

    void setLinearSolver() const override;

    void setLinearSolverB(int which) const override;

    void setNonLinearSolver() const override;

    void setNonLinearSolverSens() const override;

    void setNonLinearSolverB(int which) const override;

  protected:

    void calcIC(realtype tout1) const override;

    void calcICB(int which, realtype tout1) const override;

    void getB(int which) const override;

    void getSens() const override;

    void getQuadB(int which) const override;

    void getQuad(realtype &t_) const override;

    void getQuadDky(realtype t_, int k) const override;

    void reInitPostProcessF(realtype tnext) const override;

    void reInitPostProcessB(realtype tnext) const override;

    void reInitPostProcess(void *ami_mem, realtype *t_, AmiVector *yout,
                           realtype tout) const;

    void allocateSolver() const override;

    void setSStolerances(double rtol_, double atol_) const override;

    void setSensSStolerances(double rtol_,
                             const double *atol_) const override;

    void setSensErrCon(bool error_corr) const override;

    void setQuadErrConB(int which, bool flag) const override;

    void setQuadErrCon(bool flag) const override;

    void setErrHandlerFn() const override;

    void setUserData(Model *model) const override;

    void setUserDataB(int which, Model *model) const override;

    void setMaxNumSteps(long int mxsteps) const override;

    void setStabLimDet(int stldet_) const override;

    void setStabLimDetB(int which, int stldet_) const override;

    void setId(const Model *model) const override;

    void setSuppressAlg(bool flag) const override;

    void resetState(void *cv_mem, const_N_Vector y0) const;

    void setSensParams(const realtype *p, const realtype *pbar,
                       const int *plist) const override;

    void adjInit() const override;

    void quadInit(const AmiVector &xQ0) const override;
    
    void allocateSolverB(int *which) const override;

    void setSStolerancesB(int which, realtype relTolB,
                          realtype absTolB) const override;

    void quadSStolerancesB(int which, realtype reltolQB,
                           realtype abstolQB) const override;

    void quadSStolerances(realtype reltolQ, realtype abstolQ) const override;

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

    void getLastOrder(const void *ami_ami_mem, int *order_) const override;

    void *getAdjBmem(void *ami_mem, int which) const override;

    template <class Archive>
    friend void boost::serialization::serialize(Archive &ar, CVodeSolver &r,
                                                unsigned int version);

    friend bool operator==(const CVodeSolver &a, const CVodeSolver &b);

    void init(realtype t0, const AmiVector &x0,
              const AmiVector &dx0) const override;

    void initSteadystate(const realtype t0, const AmiVector &x0,
                         const AmiVector &dx0) const override;

    void sensInit1(const AmiVectorArray &sx0, const AmiVectorArray &sdx0)
    const override;

    void binit(int which, realtype tf, const AmiVector &xB0,
               const AmiVector &dxB0) const override;

    void qbinit(int which, const AmiVector &xQB0) const override;

    void rootInit(int ne) const override;

    void setDenseJacFn() const override;

    void setSparseJacFn() const override;

    void setBandJacFn() const override;

    void setJacTimesVecFn() const override;

    void setDenseJacFnB(int which) const override;

    void setSparseJacFnB(int which) const override;

    void setBandJacFnB(int which) const override;

    void setJacTimesVecFnB(int which) const override;

    void setSparseJacFn_ss() const override;
};

} // namespace amici

#endif /* CVodewrap_h */
