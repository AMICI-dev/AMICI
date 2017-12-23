#ifndef AMICI_SOLVER_H
#define AMICI_SOLVER_H

#include <include/amici_vector.h>
#include <include/amici_defines.h>
#include <include/amici_misc.h>
#include <include/symbolic_functions.h>
#include <nvector/nvector_serial.h>   // DlsMat
#include <sundials/sundials_sparse.h> // SlsMat

namespace amici {

class ReturnData;
class UserData;
class ForwardProblem;
class BackwardProblem;
class Model;
extern msgIdAndTxtFp warnMsgIdAndTxt;

/** Solver class.
 * provides a generic interface to CVode and IDA solvers, individual realizations
 * are realized in the CVodeSolver and the IDASolver class.
 */
class Solver {
  public:
    Solver() = default;

    virtual ~Solver() = default;

    void setupAMI(ForwardProblem *fwd, const UserData *udata, Model *model);

    void setupAMIB(BackwardProblem *bwd, const UserData *udata, Model *model);

    /**
     * AMIGetSens extracts diagnosis information from solver memory block and
     * writes them into the return data instance
     *
     * @param tret time at which the sensitivities should be computed
     * @param yySout vector with sensitivities
     */
    virtual void AMIGetSens(realtype *tret, AmiVectorArray *yySout) = 0;

    void getDiagnosis(const int it, ReturnData *rdata);

    void getDiagnosisB(const int it, ReturnData *rdata, const BackwardProblem *bwd);

    /**
     * AMIGetRootInfo extracts information which event occured
     *
     * @param rootsfound array with flags indicating whether the respective
     * event occured
     */
    virtual void AMIGetRootInfo(int *rootsfound) = 0;

    /**
     * AMIReInit reinitializes the states in the solver after an event occurence
     *
     * @param t0 new timepoint
     * @param yy0 new state variables
     * @param yp0 new derivative state variables (DAE only)
     */
    virtual void AMIReInit(realtype t0, AmiVector *yy0, AmiVector *yp0) = 0;

    /**
     * AMISensReInit reinitializes the state sensitivites in the solver after an
     * event occurence
     *
     * @param ism sensitivity mode
     * @param yS0 new state sensitivity
     * @param ypS0 new derivative state sensitivities (DAE only)
     */
    virtual void AMISensReInit(int ism, AmiVectorArray *yS0, AmiVectorArray *ypS0) = 0;

    /**
     * AMICalcIC calculates consistent initial conditions, assumes initial
     * states to be correct (DAE only)
     *
     * @param tout1 next timepoint to be computed (sets timescale)
     * @param x initial state variables
     * @param dx initial derivative state variables (DAE only)
     */
    virtual void AMICalcIC(realtype tout1, AmiVector *x, AmiVector *dx) = 0;

    /**
      * AMICalcIBC calculates consistent initial conditions for the backwards
     * problem, assumes initial states to be correct (DAE only)
      *
      * @param which identifier of the backwards problem
      * @param tout1 next timepoint to be computed (sets timescale)
      * @param xB states of final solution of the forward problem
      * @param dxB derivative states of final solution of the forward
     * problem (DAE only)
      */
    virtual void AMICalcICB(int which, realtype tout1, AmiVector *xB,
                           AmiVector *dxB) = 0;

    /**
      * AMISolve solves the forward problem until a predefined timepoint
      *
      * @param tout timepoint until which simulation should be performed
     *
      * @param yret states
      * @param ypret derivative states (DAE only)
      * @param tret pointer to the time variable
      * @param itask task identifier, can be CV_NORMAL or CV_ONE_STEP
     * @return status flag indicating success of execution
      */
    virtual int AMISolve(realtype tout, AmiVector *yret, AmiVector *ypret,
                         realtype *tret, int itask) = 0;

    /**
      * AMISolveF solves the forward problem until a predefined timepoint
     * (adjoint only)
      *
      * @param tout timepoint until which simulation should be performed
     *
      * @param yret states
      * @param ypret derivative states (DAE only)
      * @param tret pointer to the time variable
      * @param itask task identifier, can be CV_NORMAL or CV_ONE_STEP
      * @param ncheckPtr pointer to a number that counts the internal
     * checkpoints
     * @return status flag indicating success of execution
      */
    virtual int AMISolveF(realtype tout, AmiVector *yret, AmiVector *ypret,
                          realtype *tret, int itask, int *ncheckPtr) = 0;

    /**
      * AMISolveB solves the backward problem until a predefined timepoint
     * (adjoint only)
      *
      * @param tBout timepoint until which simulation should be performed
     *
      * @param itaskB task identifier, can be CV_NORMAL or CV_ONE_STEP
      */
    virtual void AMISolveB(realtype tBout, int itaskB) = 0;

    /**
      * AMISetStopTime sets a timepoint at which the simulation will be stopped
      *
      * @param tstop timepoint until which simulation should be performed
     *
      */
    virtual void AMISetStopTime(realtype tstop) = 0;

    //    virtual void AMIRootInit(int nrtfn, RootFn ptr) = 0;

    /**
      * AMIReInitB reinitializes the adjoint states after an event occurence
      *
      * @param which identifier of the backwards problem
      * @param tB0 new timepoint
      * @param yyB0 new adjoint state variables
      * @param ypB0 new adjoint derivative state variables (DAE only)
      */
    virtual void AMIReInitB(int which, realtype tB0, AmiVector *yyB0,
                           AmiVector *ypB0) = 0;

    /**
      * AMIGetB returns the current adjoint states
      *
      * @param which identifier of the backwards problem
      * @param tret time at which the adjoint states should be computed
      * @param yy adjoint state variables
      * @param yp adjoint derivative state variables (DAE only)
      */
    virtual void AMIGetB(int which, realtype *tret, AmiVector *yy,
                        AmiVector *yp) = 0;

    /**
      * AMIGetQuadB returns the current adjoint states
      *
      * @param which identifier of the backwards problem
      * @param tret time at which the adjoint states should be computed
      * @param qB adjoint quadrature state variables
      */
    virtual void AMIGetQuadB(int which, realtype *tret, AmiVector *qB) = 0;

    /**
      * AMIReInitB reinitializes the adjoint states after an event occurence
      *
      * @param which identifier of the backwards problem
      * @param yQB0 new adjoint quadrature state variables
      */
    virtual void AMIQuadReInitB(int which, AmiVector *yQB0) = 0;

    /**
      * turnOffRootFinding disables rootfinding
      */
    virtual void turnOffRootFinding() = 0;

  protected:
    /**
     * init initialises the states at the specified initial timepoint
     *
     * @param x initial state variables
     * @param dx initial derivative state variables (DAE only)
     * @param t initial timepoint
     */
    virtual void init(AmiVector *x, AmiVector *dx, realtype t) = 0;

    // TODO: check if model has adjoint sensitivities, else return -1

    /**
     * binit initialises the adjoint states at the specified final timepoint
     *
     * @param which identifier of the backwards problem
     * @param xB initial adjoint state variables
     * @param dxB initial adjoint derivative state variables (DAE only)
     * @param t final timepoint
     */
    virtual void binit(int which, AmiVector *xB, AmiVector *dxB, realtype t) = 0;

    // TODO: check if model has adjoint sensitivities, else return -1

    /**
     * qbinit initialises the quadrature states at the specified final timepoint
     *
     * @param which identifier of the backwards problem
     * @param qBdot initial adjoint quadrature state variables
     */
    virtual void qbinit(int which, AmiVector *qBdot) = 0;

    /**
     * RootInit initialises the rootfinding for events
     *
     * @param ne number of different events
     */
    virtual void rootInit(int ne) = 0;

    // TODO: check if model has forward sensitivities, else return -1
    /**
     * SensInit1 initialises the sensitivities at the specified initial
     * timepoint
     *
     * @param sx initial state sensitivities
     * @param sdx initial derivative state sensitivities (DAE only)
     * @param udata initial derivative state sensitivities (DAE only)
     */
    virtual void sensInit1(AmiVectorArray *sx, AmiVectorArray *sdx,
                          const UserData *udata) = 0;

    /**
     * SetDenseJacFn sets the dense Jacobian function
     *
     */
    virtual void setDenseJacFn() = 0;

    /**
     * SetSparseJacFn sets the sparse Jacobian function
     *
     */
    virtual void setSparseJacFn() = 0;

    /**
     * SetBandJacFn sets the banded Jacobian function
     *
     */
    virtual void setBandJacFn() = 0;

    /**
     * SetJacTimesVecFn sets the Jacobian vector multiplication function
     *
     */
    virtual void setJacTimesVecFn() = 0;

    // TODO: check if model has adjoint sensitivities, else return -1
    /**
     * SetDenseJacFn sets the dense Jacobian function
     *
     * @param which identifier of the backwards problem
     */
    virtual void setDenseJacFnB(int which) = 0;

    // TODO: check if model has adjoint sensitivities, else return -1
    /**
     * SetSparseJacFn sets the sparse Jacobian function
     *
     * @param which identifier of the backwards problem
     */
    virtual void setSparseJacFnB(int which) = 0;

    // TODO: check if model has adjoint sensitivities, else return -1
    /**
     * SetBandJacFn sets the banded Jacobian function
     *
     * @param which identifier of the backwards problem
     */
    virtual void setBandJacFnB(int which) = 0;

    // TODO: check if model has adjoint sensitivities, else return -1
    /**
     * SetJacTimesVecFn sets the Jacobian vector multiplication function
     *
     * @param which identifier of the backwards problem
     */
    virtual void setJacTimesVecFnB(int which) = 0;

    static void wrapErrHandlerFn(int error_code, const char *module,
                                 const char *function, char *msg,
                                 void *eh_data);

    /**
     * AMICreate specifies solver method and initializes solver memory for the
     * forward problem
     *
     * @param lmm linear multistep method CV_ADAMS or CV_BDF
     * @param iter nonlinear solver method CV_NEWTON or CV_FUNCTIONAL
     */
    virtual void *AMICreate(int lmm, int iter) = 0;

    /**
     * AMISStolerances sets relative and absolute tolerances for the forward
     * problem
     *
     * @param rtol relative tolerances
     * @param atol absolute tolerances
     */
    virtual void AMISStolerances(double rtol, double atol) = 0;

    /**
     * AMISensEEtolerances activates automatic estimation of tolerances for the
     * forward problem
     *
     */
    virtual void AMISensEEtolerances() = 0;

    /**
     * AMISetSensErrCon specifies whether error control is also enforced for
     * sensitivities for the forward problem
     *
     * @param error_corr activation flag
     */
    virtual void AMISetSensErrCon(bool error_corr) = 0;

    /**
     * AMISetSensErrCon specifies whether error control is also enforced for the
     * backward quadrature problem
     *
     * @param which identifier of the backwards problem
     * @param flag activation flag
     */
    virtual void AMISetQuadErrConB(int which, bool flag) = 0;

    /**
     * AMISetErrHandlerFn attaches the error handler function (errMsgIdAndTxt)
     * to the solver
     *
     */
    virtual void AMISetErrHandlerFn() = 0;

    /**
     * AMISetUserData attaches the user data instance (here this is a TempData and
     * not UserData instance) to the forward problem
     *
     * @param model Model instance,
     */
    virtual void AMISetUserData(Model *model) = 0;

    /**
     * AMISetUserDataB attaches the user data instance (here this is a TempData
     * and not UserData instance) to the backward problem
     *
     * @param which identifier of the backwards problem
     * @param model Model instance,
     */
    virtual void AMISetUserDataB(int which, Model *model) = 0;

    /**
     * AMISetMaxNumSteps specifies the maximum number of steps for the forward
     * problem
     *
     * @param mxsteps number of steps
     */
    virtual void AMISetMaxNumSteps(long int mxsteps) = 0;

    /**
     * AMISetMaxNumStepsB specifies the maximum number of steps for the forward
     * problem
     *
     * @param which identifier of the backwards problem
     * @param mxstepsB number of steps
     */
    virtual void AMISetMaxNumStepsB(int which, long int mxstepsB) = 0;

    /**
     * AMISetStabLimDet activates stability limit detection for the forward
     * problem
     *
     * @param stldet flag for stability limit detection (TRUE or FALSE)
     *
     */
    virtual void AMISetStabLimDet(int stldet) = 0;

    /**
     * AMISetStabLimDetB activates stability limit detection for the backward
     * problem
     *
     * @param which identifier of the backwards problem
     * @param stldet flag for stability limit detection (TRUE or FALSE)
     *
     */
    virtual void AMISetStabLimDetB(int which, int stldet) = 0;

    /**
     * AMISetId specify algebraic/differential components (DAE only)
     *
     * @param model model specification
     */
    virtual void AMISetId(Model *model) = 0;

    /**
     * AMISetId deactivates error control for algebraic components (DAE only)
     *
     * @param flag deactivation flag
     */
    virtual void AMISetSuppressAlg(bool flag) = 0;

    /**
     * AMISetSensParams specifies the scaling and indexes for sensitivity
     * computation
     *
     * @param p paramaters
     * @param pbar parameter scaling constants
     * @param plist parameter index list
     */
    virtual void AMISetSensParams(realtype *p, realtype *pbar, int *plist) = 0;

    /**
     * AMIGetDky interpolates the (derivative of the) solution at the requested
     * timepoint
     *
     * @param t timepoint
     * @param k derivative order
     * @param dky interpolated solution
     */
    virtual void AMIGetDky(realtype t, int k, AmiVector *dky) = 0;

    /**
     * AMIFree frees allocation solver memory
     */
    virtual void AMIFree() = 0;

    /**
     * AMIAdjInit initializes the adjoint problem
     *
     * @param steps number of integration points between checkpoints
     * @param interp interpolation type, can be CV_POLYNOMIAL or CV_HERMITE
     *
     */
    virtual void AMIAdjInit(long int steps, int interp) = 0;

    /**
     * AMICreateB specifies solver method and initializes solver memory for the
     * backward problem
     *
     * @param which identifier of the backwards problem
     * @param lmm linear multistep method CV_ADAMS or CV_BDF
     * @param iter nonlinear solver method CV_NEWTON or CV_FUNCTIONAL
     */
    virtual void AMICreateB(int lmm, int iter, int *which) = 0;

    /**
     * AMISStolerancesB sets relative and absolute tolerances for the backward
     * problem
     *
     * @param which identifier of the backwards problem
     * @param relTolB relative tolerances
     * @param absTolB absolute tolerances
     */
    virtual void AMISStolerancesB(int which, realtype relTolB,
                                 realtype absTolB) = 0;

    /**
     * AMISStolerancesB sets relative and absolute tolerances for the quadrature
     * backward problem
     *
     * @param which identifier of the backwards problem
     * @param reltolQB relative tolerances
     * @param abstolQB absolute tolerances
     */
    virtual void AMIQuadSStolerancesB(int which, realtype reltolQB,
                                     realtype abstolQB) = 0;

    /**
     * AMIDense attaches a dense linear solver to the forward problem
     *
     * @param nx number of state variables
     */
    virtual void AMIDense(int nx) = 0;

    /**
     * AMIDenseB attaches a dense linear solver to the backward problem
     *
     * @param which identifier of the backwards problem
     * @param nx number of state variables
     */
    virtual void AMIDenseB(int which, int nx) = 0;

    /**
     * AMIBand attaches a banded linear solver to the forward problem
     *
     * @param nx number of state variables
     * @param ubw upper matrix bandwidth
     * @param lbw lower matrix bandwidth
     */
    virtual void AMIBand(int nx, int ubw, int lbw) = 0;

    /**
     * AMIBandB attaches a banded linear solver to the backward problem
     *
     * @param which identifier of the backwards problem
     * @param nx number of state variables
     * @param ubw upper matrix bandwidth
     * @param lbw lower matrix bandwidth
     */
    virtual void AMIBandB(int which, int nx, int ubw, int lbw) = 0;

    /**
     * AMIDiag attaches a diagonal linear solver to the forward problem
     *
     */
    virtual void AMIDiag() = 0;

    /**
     * AMIDiagB attaches a diagonal linear solver to the backward problem
     *
     * @param which identifier of the backwards problem
     */
    virtual void AMIDiagB(int which) = 0;

    /**
     * AMIDAMISpgmr attaches a scaled predonditioned GMRES linear solver to the
     * forward problem
     *
     * @param prectype preconditioner type PREC_NONE, PREC_LEFT, PREC_RIGHT
     * or PREC_BOTH
     * @param maxl maximum Kryloc subspace dimension
     */
    virtual void AMISpgmr(int prectype, int maxl) = 0;

    /**
     * AMIDAMISpgmrB attaches a scaled predonditioned GMRES linear solver to the
     * backward problem
     *
     * @param which identifier of the backwards problem
     * @param prectype preconditioner type PREC_NONE, PREC_LEFT, PREC_RIGHT
     * or PREC_BOTH
     * @param maxl maximum Kryloc subspace dimension
     */
    virtual void AMISpgmrB(int which, int prectype, int maxl) = 0;

    /**
     * AMISpbcg attaches a scaled predonditioned Bi-CGStab linear solver to the
     * forward problem
     *
     * @param prectype preconditioner type PREC_NONE, PREC_LEFT, PREC_RIGHT
     * or PREC_BOTH
     * @param maxl maximum Kryloc subspace dimension
     */
    virtual void AMISpbcg(int prectype, int maxl) = 0;

    /**
     * AMISpbcgB attaches a scaled predonditioned Bi-CGStab linear solver to the
     * backward problem
     *
     * @param which identifier of the backwards problem
     * @param prectype preconditioner type PREC_NONE, PREC_LEFT, PREC_RIGHT
     * or PREC_BOTH
     * @param maxl maximum Kryloc subspace dimension
     */
    virtual void AMISpbcgB(int which, int prectype, int maxl) = 0;

    /**
     * AMISptfqmr attaches a scaled predonditioned TFQMR linear solver to the
     * forward problem
     *
     * @param prectype preconditioner type PREC_NONE, PREC_LEFT, PREC_RIGHT
     * or PREC_BOTH
     * @param maxl maximum Kryloc subspace dimension
     */
    virtual void AMISptfqmr(int prectype, int maxl) = 0;

    /**
     * AMISptfqmrB attaches a scaled predonditioned TFQMR linear solver to the
     * backward problem
     *
     * @param which identifier of the backwards problem
     * @param prectype preconditioner type PREC_NONE, PREC_LEFT, PREC_RIGHT
     * or PREC_BOTH
     * @param maxl maximum Kryloc subspace dimension
     */
    virtual void AMISptfqmrB(int which, int prectype, int maxl) = 0;

    /**
     * AMIKLU attaches a sparse linear solver to the forward problem
     *
     * @param nx number of state variables
     * @param nnz number of nonzero entries in the jacobian
     * @param sparsetype sparse storage type, CSC_MAT for column matrix,
     * CSR_MAT for row matrix
     */
    virtual void AMIKLU(int nx, int nnz, int sparsetype) = 0;

    /**
     * AMIKLUSetOrdering sets the ordering for the sparse linear solver of the
     * forward problem
     *
     * @param ordering ordering algorithm to reduce fill 0:AMD 1:COLAMD 2:
     * natural ordering
     */
    virtual void AMIKLUSetOrdering(int ordering) = 0;

    /**
     * AMIKLUSetOrderingB sets the ordering for the sparse linear solver of the
     * backward problem
     *
     * @param which identifier of the backwards problem
     * @param ordering ordering algorithm to reduce fill 0:AMD 1:COLAMD 2:
     * natural ordering
     */
    virtual void AMIKLUSetOrderingB(int which, int ordering) = 0;

    /**
     * AMIKLUB attaches a sparse linear solver to the forward problem
     *
     * @param which identifier of the backwards problem
     * @param nx number of state variables
     * @param nnz number of nonzero entries in the jacobian
     * @param sparsetype sparse storage type, CSC_MAT for column matrix,
     * CSR_MAT for row matrix
     */
    virtual void AMIKLUB(int which, int nx, int nnz, int sparsetype) = 0;

    /**
     * AMIGetNumSteps reports the number of solver steps
     *
     * @param ami_mem pointer to the solver memory instance (can be from
     * forward or backward problem)
     * @param numsteps output array
     */
    virtual void AMIGetNumSteps(void *ami_mem, long int *numsteps) = 0;

    /**
     * AMIGetNumRhsEvals reports the number of right hand evaluations
     *
     * @param ami_mem pointer to the solver memory instance (can be from
     * forward or backward problem)
     * @param numrhsevals output array
     */
    virtual void AMIGetNumRhsEvals(void *ami_mem, long int *numrhsevals) = 0;

    /**
     * AMIGetNumErrTestFails reports the number of local error test failures
     *
     * @param ami_mem pointer to the solver memory instance (can be from
     * forward or backward problem)
     * @param numerrtestfails output array
     */
    virtual void AMIGetNumErrTestFails(void *ami_mem,
                                      long int *numerrtestfails) = 0;

    /**
     * AMIGetNumNonlinSolvConvFails reports the number of nonlinear convergence
     * failures
     *
     * @param ami_mem pointer to the solver memory instance (can be from
     * forward or backward problem)
     * @param numnonlinsolvconvfails output array
     */
    virtual void
    AMIGetNumNonlinSolvConvFails(void *ami_mem,
                                 long int *numnonlinsolvconvfails) = 0;

    /**
     * AMIGetLastOrder reports the order of the integration method during the
     * last internal step
     *
     * @param ami_mem pointer to the solver memory instance (can be from
     * forward or backward problem)
     * @param order output array
     */
    virtual void AMIGetLastOrder(void *ami_mem, int *order) = 0;

    /**
     * AMIGetAdjBmem retrieves the solver memory instance for the backward problem
     *
     * @param which identifier of the backwards problem
     * @param ami_mem pointer to the forward solver memory instance
     * @return ami_memB pointer to the backward solver memory instance
     */
    virtual void *AMIGetAdjBmem(void *ami_mem, int which) = 0;

    void setLinearSolver(const UserData *udata, Model *model);
    void setLinearSolverB(const UserData *udata, Model *model, int which);

protected:
    
    /** pointer to ami memory block */
    void *ami_mem = nullptr;
    
    bool solverWasCalled = false;
};

} // namespace amici

#endif // AMICISOLVER_H
