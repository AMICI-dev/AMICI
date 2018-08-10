#ifndef AMICI_SOLVER_CVODES_h
#define AMICI_SOLVER_CVODES_h

#include "amici/solver.h"
#include "amici/defines.h"
#include "amici/vector.h"

#include <cvodes/cvodes_dense.h>
#include <sundials/sundials_sparse.h>

namespace amici {

class ExpData;
class ReturnData;
class Model_ODE;
class CVodeSolver;
}

// for serialization friend in Solver
namespace boost { namespace serialization {
template <class Archive>
void serialize(Archive &ar, amici::CVodeSolver &u, const unsigned int version);
}}


namespace amici {

class CVodeSolver : public Solver {
  public:
    CVodeSolver() = default;

    ~CVodeSolver() override = default;

    /**
     * @brief Clone this instance
     * @return The clone
     */
    virtual Solver* clone() const override;
    
    void reInit(realtype t0, AmiVector *yy0, AmiVector *yp0) override;
    
    void sensReInit( AmiVectorArray *yS0, AmiVectorArray *ypS0) override;
    
    void reInitB(int which, realtype tB0, AmiVector *yyB0,
                 AmiVector *ypB0) override;
    
    void quadReInitB(int which, AmiVector *yQB0) override;
    
    int solve(realtype tout, AmiVector *yret, AmiVector *ypret, realtype *tret,
              int itask) override;
    
    int solveF(realtype tout, AmiVector *yret, AmiVector *ypret, realtype *tret,
               int itask, int *ncheckPtr) override;
    
    void solveB(realtype tBout, int itaskB) override;
    
    void getB(int which, realtype *tret, AmiVector *yy, AmiVector *yp) const override;
    
    void getDky(realtype t, int k, AmiVector *dky) const override;
    
    void getSens(realtype *tret, AmiVectorArray *yySout) const override;
    
    void getQuadB(int which, realtype *tret, AmiVector *qB) const override;
    
    void getRootInfo(int *rootsfound) const override;
    
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
    
    static int fxdot(realtype t, N_Vector x, N_Vector xdot, void *user_data);
    
    static int fJSparse(realtype t, N_Vector x, N_Vector xdot, SlsMat J,
                        void *user_data, N_Vector tmp1, N_Vector tmp2,
                        N_Vector tmp3);
    
    static int fJ(long int N, realtype t, N_Vector x, N_Vector xdot,
                  DlsMat J, void *user_data, N_Vector tmp1,
                  N_Vector tmp2, N_Vector tmp3);
    
protected:
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

    void setSensParams(realtype *p, realtype *pbar, int *plist) override;

    void adjInit() override;

    void allocateSolverB(int *which) override;

    void setSStolerancesB(int which, realtype relTolB,
                         realtype absTolB) override;

    void quadSStolerancesB(int which, realtype reltolQB,
                             realtype abstolQB) override;

    void setMaxNumStepsB(int which, long int mxstepsB) override;

    void dense(int nx) override;

    void denseB(int which, int nx) override;

    void band(int nx, int ubw, int lbw) override;

    void bandB(int which, int nx, int ubw, int lbw) override;

    void diag() override;

    void diagB(int which) override;

    void spgmr(int prectype, int maxl) override;

    void spgmrB(int which, int prectype, int maxl) override;

    void spbcg(int prectype, int maxl) override;

    void spbcgB(int which, int prectype, int maxl) override;

    void sptfqmr(int prectype, int maxl) override;

    void sptfqmrB(int which, int prectype, int maxl) override;

    void klu(int nx, int nnz, int sparsetype) override;

    void kluSetOrdering(int ordering) override;

    void kluSetOrderingB(int which, int ordering) override;

    void kluB(int which, int nx, int nnz, int sparsetype) override;

    void getNumSteps(void *ami_mem, long int *numsteps) const override;

    void getNumRhsEvals(void *ami_mem, long int *numrhsevals) const override;

    void getNumErrTestFails(void *ami_mem,
                              long int *numerrtestfails) const override;

    void getNumNonlinSolvConvFails(void *ami_mem,
                                     long int *numnonlinsolvconvfails) const override;

    void getLastOrder(void *ami_ami_mem, int *order) const override;

    void *getAdjBmem(void *ami_mem, int which) override;
    
    template <class Archive>
    friend void boost::serialization::serialize(Archive &ar, CVodeSolver &r, const unsigned int version);

    friend bool operator ==(const CVodeSolver &a, const CVodeSolver &b);
    
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
    
    static int fJB(long int NeqBdot, realtype t, N_Vector x, N_Vector xB,
                   N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B,
                   N_Vector tmp2B, N_Vector tmp3B);
    
    static int fJSparseB(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot,
                         SlsMat JB, void *user_data, N_Vector tmp1B,
                         N_Vector tmp2B, N_Vector tmp3B);
    
    static int fJBand(long int N, long int mupper, long int mlower, realtype t,
                      N_Vector x, N_Vector xdot, DlsMat J, void *user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    
    static int fJBandB(long int NeqBdot, long int mupper, long int mlower,
                       realtype t, N_Vector x, N_Vector xB, N_Vector xBdot,
                       DlsMat JB, void *user_data, N_Vector tmp1B,
                       N_Vector tmp2B, N_Vector tmp3B);
    
    static int fJDiag(realtype t, N_Vector JDiag, N_Vector x, void *user_data);
    
    static int fJv(N_Vector v, N_Vector Jv, realtype t, N_Vector x, N_Vector xdot,
                   void *user_data, N_Vector tmp);
    
    static int fJvB(N_Vector vB, N_Vector JvB, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot,
                    void *user_data, N_Vector tmpB);
    
    static int froot(realtype t, N_Vector x, realtype *root,
                     void *user_data);
    
    static int fxBdot(realtype t, N_Vector x, N_Vector xB,
                      N_Vector xBdot, void *user_data);
    
    static int fqBdot(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot,
                      void *user_data);
    
    static int fsxdot(int Ns, realtype t, N_Vector x, N_Vector xdot, int ip,
                      N_Vector sx, N_Vector sxdot, void *user_data,
                      N_Vector tmp1, N_Vector tmp2);
};

} // namespace amici

#endif /* CVodewrap_h */
