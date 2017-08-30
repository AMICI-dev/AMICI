#ifndef AMICISOLVER_H
#define AMICISOLVER_H

#include <nvector/nvector_serial.h>   // DlsMat
#include <sundials/sundials_sparse.h> // SlsMat

class ReturnData;
class UserData;
class TempData;
class Model;

//!  Solver class.
/*!
  provides a generic interface to CVode and IDA solvers, individual realizations
  are realized in the CVodeSolver and the IDASolver class.
*/
class Solver {
  public:
    Solver() {}

    virtual ~Solver();

    int setupAMI(const UserData *udata, TempData *tdata, Model *model);

    int setupAMIB(const UserData *udata, TempData *tdata, Model *model);

    /**
     * AMIGetSens extracts diagnosis information from solver memory block and
     * writes them into the return data object
     *
     * @param[in] tret time at which the sensitivities should be computed
     * @param[out] yySout vector with sensitivities @type N_Vector
     * @return status flag indicating success of execution @type int
     */
    virtual int AMIGetSens(realtype *tret, N_Vector *yySout) = 0;

    int getDiagnosis(const int it, ReturnData *rdata);

    int getDiagnosisB(const int it, ReturnData *rdata, const TempData *tdata);

    /**
     * AMIGetRootInfo extracts information which event occured
     *
     * @param[out] rootsfound array with flags indicating whether the respective
     * event occured
     * @return status flag indicating success of execution @type int
     */
    virtual int AMIGetRootInfo(int *rootsfound) = 0;

    /**
     * AMIReInit reinitializes the states in the solver after an event occurence
     *
     * @param[in] t0 new timepoint @type realtype
     * @param[in] yy0 new state variables @type N_Vector
     * @param[in] yp0 new derivative state variables (DAE only) @type N_Vector
     * @return status flag indicating success of execution @type int
     */
    virtual int AMIReInit(realtype t0, N_Vector yy0, N_Vector yp0) = 0;

    /**
     * AMISensReInit reinitializes the state sensitivites in the solver after an
     * event occurence
     *
     * @param[in] ism sensitivity mode @type realtype
     * @param[in] yS0 new state sensitivity @type N_Vector
     * @param[in] ypS0 new derivative state sensitivities (DAE only) @type
     * N_Vector
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISensReInit(int ism, N_Vector *yS0, N_Vector *ypS0) = 0;

    /**
      * AMICalcIC calculates consistent initial conditions, assumes initial
     * states to be correct (DAE only)
      *
      * @param[in] tout1 next timepoint to be computed (sets timescale) @type
     * realtype
      * @return status flag indicating success of execution @type int
      */
    virtual int AMICalcIC(realtype tout1) = 0;

    /**
      * AMICalcIBC calculates consistent initial conditions for the backwards
     * problem, assumes initial states to be correct (DAE only)
      *
      * @param[in] which identifier of the backwards problem @type int
      * @param[in] tout1 next timepoint to be computed (sets timescale) @type
     * realtype
      * @param[in] xB states of final solution of the forward problem @type
     * N_Vector
      * @param[in] dxB derivative states of final solution of the forward
     * problem (DAE only) @type N_Vector
      * @return status flag indicating success of execution @type int
      */
    virtual int AMICalcICB(int which, realtype tout1, N_Vector xB,
                           N_Vector dxB) = 0;

    /**
      * AMISolve solves the forward problem until a predefined timepoint
      *
      * @param[in] tout timepoint until which simulation should be performed
     * @type realtype
      * @param[in] yret states @type N_Vector
      * @param[in] ypret derivative states (DAE only) @type N_Vector
      * @param[in,out] tret pointer to the time variable @type realtype
      * @param[in] itask task identifier, can be CV_NORMAL or CV_ONE_STEP @type
     * realtype
      * @return status flag indicating success of execution @type int
      */
    virtual int AMISolve(realtype tout, N_Vector yret, N_Vector ypret,
                         realtype *tret, int itask) = 0;

    /**
      * AMISolveF solves the forward problem until a predefined timepoint
     * (adjoint only)
      *
      * @param[in] tout timepoint until which simulation should be performed
     * @type realtype
      * @param[in] yret states @type N_Vector
      * @param[in] ypret derivative states (DAE only) @type N_Vector
      * @param[in,out] tret pointer to the time variable @type realtype
      * @param[in] itask task identifier, can be CV_NORMAL or CV_ONE_STEP @type
     * realtype
      * @param[in] ncheckPtr pointer to a number that counts the internal
     * checkpoints @type realtype
      * @return status flag indicating success of execution @type int
      */
    virtual int AMISolveF(realtype tout, N_Vector yret, N_Vector ypret,
                          realtype *tret, int itask, int *ncheckPtr) = 0;

    /**
      * AMISolveB solves the backward problem until a predefined timepoint
     * (adjoint only)
      *
      * @param[in] tBout timepoint until which simulation should be performed
     * @type realtype
      * @param[in] itaskB task identifier, can be CV_NORMAL or CV_ONE_STEP @type
     * realtype
      * @return status flag indicating success of execution @type int
      */
    virtual int AMISolveB(realtype tBout, int itaskB) = 0;

    /**
      * AMISetStopTime sets a timepoint at which the simulation will be stopped
      *
      * @param[in] tstop timepoint until which simulation should be performed
     * @type realtype
      * @return status flag indicating success of execution @type int
      */
    virtual int AMISetStopTime(realtype tstop) = 0;

    //    virtual int AMIRootInit(int nrtfn, RootFn ptr) = 0;

    /**
      * AMIReInitB reinitializes the adjoint states after an event occurence
      *
      * @param[in] which identifier of the backwards problem @type int
      * @param[in] tB0 new timepoint @type realtype
      * @param[in] yyB0 new adjoint state variables @type N_Vector
      * @param[in] ypB0 new adjoint derivative state variables (DAE only) @type
     * N_Vector
      * @return status flag indicating success of execution @type int
      */
    virtual int AMIReInitB(int which, realtype tB0, N_Vector yyB0,
                           N_Vector ypB0) = 0;

    /**
      * AMIGetB returns the current adjoint states
      *
      * @param[in] which identifier of the backwards problem @type int
      * @param[in] tret time at which the adjoint states should be computed
      * @param[in] yy adjoint state variables @type N_Vector
      * @param[in] yp adjoint derivative state variables (DAE only) @type
     * N_Vector
      * @return status flag indicating success of execution @type int
      */
    virtual int AMIGetB(int which, realtype *tret, N_Vector yy,
                        N_Vector yp) = 0;

    /**
      * AMIGetQuadB returns the current adjoint states
      *
      * @param[in] which identifier of the backwards problem @type int
      * @param[in] tret time at which the adjoint states should be computed
      * @param[in] qB adjoint quadrature state variables @type N_Vector
      * @return status flag indicating success of execution @type int
      */
    virtual int AMIGetQuadB(int which, realtype *tret, N_Vector qB) = 0;

    /**
      * AMIReInitB reinitializes the adjoint states after an event occurence
      *
      * @param[in] which identifier of the backwards problem @type int
      * @param[in] yQB0 new adjoint quadrature state variables @type N_Vector
      * @return status flag indicating success of execution @type int
      */
    virtual int AMIQuadReInitB(int which, N_Vector yQB0) = 0;

    /**
      * turnOffRootFinding disables rootfinding
      *
      * @return status flag indicating success of execution @type int
      */
    virtual int turnOffRootFinding() = 0;

  protected:
    /**
     * init initialises the states at the specified initial timepoint
     *
     * @param[in] x initial state variables @type N_Vector
     * @param[in] dx initial derivative state variables (DAE only) @type
     * N_Vector
     * @param[in] t initial timepoint @type realtype
     * @return status flag indicating success of execution @type int
     */
    virtual int init(N_Vector x, N_Vector dx, realtype t) = 0;

    // TODO: check if model has adjoint sensitivities, else return -1

    /**
     * binit initialises the adjoint states at the specified final timepoint
     *
     * @param[in] which identifier of the backwards problem @type int
     * @param[in] xB initial adjoint state variables @type N_Vector
     * @param[in] dxB initial adjoint derivative state variables (DAE only)
     * @type N_Vector
     * @param[in] t final timepoint @type realtype
     * @return status flag indicating success of execution @type int
     */
    virtual int binit(int which, N_Vector xB, N_Vector dxB, realtype t) = 0;

    // TODO: check if model has adjoint sensitivities, else return -1

    /**
     * qbinit initialises the quadrature states at the specified final timepoint
     *
     * @param[in] which identifier of the backwards problem @type int
     * @param[in] qBdot initial adjoint quadrature state variables @type
     * N_Vector
     * @return status flag indicating success of execution @type int
     */
    virtual int qbinit(int which, N_Vector qBdot) = 0;

    /**
     * RootInit initialises the rootfinding for events
     *
     * @param[in] ne number of different events @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int rootInit(int ne) = 0;

    // TODO: check if model has forward sensitivities, else return -1
    /**
     * SensInit1 initialises the sensitivities at the specified initial
     * timepoint
     *
     * @param[in] sx initial state sensitivities @type N_Vector
     * @param[in] sdx initial derivative state sensitivities (DAE only) @type
     * N_Vector
     * @param[in] udata initial derivative state sensitivities (DAE only) @type
     * N_Vector
     * @return status flag indicating success of execution @type int
     */
    virtual int sensInit1(N_Vector *sx, N_Vector *sdx,
                          const UserData *udata) = 0;

    /**
     * SetDenseJacFn sets the dense Jacobian function
     *
     * @return status flag indicating success of execution @type int
     */
    virtual int setDenseJacFn() = 0;

    /**
     * SetSparseJacFn sets the sparse Jacobian function
     *
     * @return status flag indicating success of execution @type int
     */
    virtual int setSparseJacFn() = 0;

    /**
     * SetBandJacFn sets the banded Jacobian function
     *
     * @return status flag indicating success of execution @type int
     */
    virtual int setBandJacFn() = 0;

    /**
     * SetJacTimesVecFn sets the Jacobian vector multiplication function
     *
     * @return status flag indicating success of execution @type int
     */
    virtual int setJacTimesVecFn() = 0;

    // TODO: check if model has adjoint sensitivities, else return -1
    /**
     * SetDenseJacFn sets the dense Jacobian function
     *
     * @param[in] which identifier of the backwards problem @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int setDenseJacFnB(int which) = 0;

    // TODO: check if model has adjoint sensitivities, else return -1
    /**
     * SetSparseJacFn sets the sparse Jacobian function
     *
     * @param[in] which identifier of the backwards problem @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int setSparseJacFnB(int which) = 0;

    // TODO: check if model has adjoint sensitivities, else return -1
    /**
     * SetBandJacFn sets the banded Jacobian function
     *
     * @param[in] which identifier of the backwards problem @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int setBandJacFnB(int which) = 0;

    // TODO: check if model has adjoint sensitivities, else return -1
    /**
     * SetJacTimesVecFn sets the Jacobian vector multiplication function
     *
     * @param[in] which identifier of the backwards problem @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int setJacTimesVecFnB(int which) = 0;

    static void wrapErrHandlerFn(int error_code, const char *module,
                                 const char *function, char *msg,
                                 void *eh_data);

    /**
     * AMICreate specifies solver method and initializes solver memory for the
     * forward problem
     *
     * @param[in] lmm linear multistep method CV_ADAMS or CV_BDF @type int
     * @param[in] iter nonlinear solver method CV_NEWTON or CV_FUNCTIONAL @type
     * int
     * @return status flag indicating success of execution @type int
     */
    virtual void *AMICreate(int lmm, int iter) = 0;

    /**
     * AMISStolerances sets relative and absolute tolerances for the forward
     * problem
     *
     * @param[in] rtol relative tolerances @type double
     * @param[in] atol absolute tolerances @type double
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISStolerances(double rtol, double atol) = 0;

    /**
     * AMISensEEtolerances activates automatic estimation of tolerances for the
     * forward problem
     *
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISensEEtolerances() = 0;

    /**
     * AMISetSensErrCon specifies whether error control is also enforced for
     * sensitivities for the forward problem
     *
     * @param[in] error_corr activation flag @type bool
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISetSensErrCon(bool error_corr) = 0;

    /**
     * AMISetSensErrCon specifies whether error control is also enforced for the
     * backward quadrature problem
     *
     * @param[in] which identifier of the backwards problem @type int
     * @param[in] flag activation flag @type bool
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISetQuadErrConB(int which, bool flag) = 0;

    /**
     * AMISetErrHandlerFn attaches the error handler function (errMsgIdAndTxt)
     * to the solver
     *
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISetErrHandlerFn() = 0;

    /**
     * AMISetUserData attaches the user data object (here this is a TempData and
     * not UserData object) to the forward problem
     *
     * @param[in] user_data TempData object, @type TempData
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISetUserData(void *user_data) = 0;

    /**
     * AMISetUserDataB attaches the user data object (here this is a TempData
     * and not UserData object) to the backward problem
     *
     * @param[in] which identifier of the backwards problem @type int
     * @param[in] user_data TempData object, @type TempData
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISetUserDataB(int which, void *user_data) = 0;

    /**
     * AMISetMaxNumSteps specifies the maximum number of steps for the forward
     * problem
     *
     * @param[in] mxsteps number of steps @type long int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISetMaxNumSteps(long int mxsteps) = 0;

    /**
     * AMISetMaxNumStepsB specifies the maximum number of steps for the forward
     * problem
     *
     * @param[in] which identifier of the backwards problem @type int
     * @param[in] mxstepsB number of steps @type long int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISetMaxNumStepsB(int which, long int mxstepsB) = 0;

    /**
     * AMISetStabLimDet activates stability limit detection for the forward
     * problem
     *
     * @param[in] stldet flag for stability limit detection (TRUE or FALSE)
     * @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISetStabLimDet(int stldet) = 0;

    /**
     * AMISetStabLimDetB activates stability limit detection for the backward
     * problem
     *
     * @param[in] which identifier of the backwards problem @type int
     * @param[in] stldet flag for stability limit detection (TRUE or FALSE)
     * @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISetStabLimDetB(int which, int stldet) = 0;

    /**
     * AMISetId specify algebraic/differential components (DAE only)
     *
     * @param[in] model model specification @type Model
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISetId(Model *model) = 0;

    /**
     * AMISetId deactivates error control for algebraic components (DAE only)
     *
     * @param[in] flag deactivation flag @type bool
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISetSuppressAlg(bool flag) = 0;

    /**
     * AMISetSensParams specifies the scaling and indexes for sensitivity
     * computation
     *
     * @param[in] p paramaters @type realtype
     * @param[in] pbar parameter scaling constants @type realtype
     * @param[in] plist parameter index list @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISetSensParams(realtype *p, realtype *pbar, int *plist) = 0;

    /**
     * AMIGetDky interpolates the (derivative of the) solution at the requested
     * timepoint
     *
     * @param[in] t timepoint @type realtype
     * @param[in] k derivative order @type int
     * @param[out] dky interpolated solution @type N_Vector
     * @return status flag indicating success of execution @type int
     */
    virtual int AMIGetDky(realtype t, int k, N_Vector dky) = 0;

    /**
     * AMIFree frees allocation solver memory
     *
     * @return status flag indicating success of execution @type int
     */
    virtual void AMIFree() = 0;

    /**
     * AMIAdjInit initializes the adjoint problem
     *
     * @param[in] steps number of integration points between checkpoints @type
     * long int
     * @param[in] interp interpolation type, can be CV_POLYNOMIAL or CV_HERMITE
     * @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMIAdjInit(long int steps, int interp) = 0;

    /**
     * AMICreateB specifies solver method and initializes solver memory for the
     * backward problem
     *
     * @param[in] which identifier of the backwards problem @type int
     * @param[in] lmm linear multistep method CV_ADAMS or CV_BDF @type int
     * @param[in] iter nonlinear solver method CV_NEWTON or CV_FUNCTIONAL @type
     * int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMICreateB(int lmm, int iter, int *which) = 0;

    /**
     * AMISStolerancesB sets relative and absolute tolerances for the backward
     * problem
     *
     * @param[in] which identifier of the backwards problem @type int
     * @param[in] relTolB relative tolerances @type double
     * @param[in] absTolB absolute tolerances @type double
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISStolerancesB(int which, realtype relTolB,
                                 realtype absTolB) = 0;

    /**
     * AMISStolerancesB sets relative and absolute tolerances for the quadrature
     * backward problem
     *
     * @param[in] which identifier of the backwards problem @type int
     * @param[in] reltolQB relative tolerances @type double
     * @param[in] abstolQB absolute tolerances @type double
     * @return status flag indicating success of execution @type int
     */
    virtual int AMIQuadSStolerancesB(int which, realtype reltolQB,
                                     realtype abstolQB) = 0;

    /**
     * AMIDense attaches a dense linear solver to the forward problem
     *
     * @param[in] nx number of state variables @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMIDense(int nx) = 0;

    /**
     * AMIDenseB attaches a dense linear solver to the backward problem
     *
     * @param[in] which identifier of the backwards problem @type int
     * @param[in] nx number of state variables @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMIDenseB(int which, int nx) = 0;

    /**
     * AMIBand attaches a banded linear solver to the forward problem
     *
     * @param[in] nx number of state variables @type int
     * @param[in] ubw upper matrix bandwidth @type int
     * @param[in] lbw lower matrix bandwidth @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMIBand(int nx, int ubw, int lbw) = 0;

    /**
     * AMIBandB attaches a banded linear solver to the backward problem
     *
     * @param[in] which identifier of the backwards problem @type int
     * @param[in] nx number of state variables @type int
     * @param[in] ubw upper matrix bandwidth @type int
     * @param[in] lbw lower matrix bandwidth @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMIBandB(int which, int nx, int ubw, int lbw) = 0;

    /**
     * AMIDiag attaches a diagonal linear solver to the forward problem
     *
     * @return status flag indicating success of execution @type int
     */
    virtual int AMIDiag() = 0;

    /**
     * AMIDiagB attaches a diagonal linear solver to the backward problem
     *
     * @param[in] which identifier of the backwards problem @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMIDiagB(int which) = 0;

    /**
     * AMIDAMISpgmr attaches a scaled predonditioned GMRES linear solver to the
     * forward problem
     *
     * @param[in] prectype preconditioner type PREC_NONE, PREC_LEFT, PREC_RIGHT
     * or PREC_BOTH @type int
     * @param[in] maxl maximum Kryloc subspace dimension @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISpgmr(int prectype, int maxl) = 0;

    /**
     * AMIDAMISpgmrB attaches a scaled predonditioned GMRES linear solver to the
     * backward problem
     *
     * @param[in] which identifier of the backwards problem @type int
     * @param[in] prectype preconditioner type PREC_NONE, PREC_LEFT, PREC_RIGHT
     * or PREC_BOTH @type int
     * @param[in] maxl maximum Kryloc subspace dimension @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISpgmrB(int which, int prectype, int maxl) = 0;

    /**
     * AMISpbcg attaches a scaled predonditioned Bi-CGStab linear solver to the
     * forward problem
     *
     * @param[in] prectype preconditioner type PREC_NONE, PREC_LEFT, PREC_RIGHT
     * or PREC_BOTH @type int
     * @param[in] maxl maximum Kryloc subspace dimension @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISpbcg(int prectype, int maxl) = 0;

    /**
     * AMISpbcgB attaches a scaled predonditioned Bi-CGStab linear solver to the
     * backward problem
     *
     * @param[in] which identifier of the backwards problem @type int
     * @param[in] prectype preconditioner type PREC_NONE, PREC_LEFT, PREC_RIGHT
     * or PREC_BOTH @type int
     * @param[in] maxl maximum Kryloc subspace dimension @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISpbcgB(int which, int prectype, int maxl) = 0;

    /**
     * AMISptfqmr attaches a scaled predonditioned TFQMR linear solver to the
     * forward problem
     *
     * @param[in] prectype preconditioner type PREC_NONE, PREC_LEFT, PREC_RIGHT
     * or PREC_BOTH @type int
     * @param[in] maxl maximum Kryloc subspace dimension @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISptfqmr(int prectype, int maxl) = 0;

    /**
     * AMISptfqmrB attaches a scaled predonditioned TFQMR linear solver to the
     * backward problem
     *
     * @param[in] which identifier of the backwards problem @type int
     * @param[in] prectype preconditioner type PREC_NONE, PREC_LEFT, PREC_RIGHT
     * or PREC_BOTH @type int
     * @param[in] maxl maximum Kryloc subspace dimension @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISptfqmrB(int which, int prectype, int maxl) = 0;

    /**
     * AMIKLU attaches a sparse linear solver to the forward problem
     *
     * @param[in] nx number of state variables @type int
     * @param[in] nnz number of nonzero entries in the jacobian @type int
     * @param[in] sparsetype sparse storage type, CSC_MAT for column matrix,
     * CSR_MAT for row matrix @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMIKLU(int nx, int nnz, int sparsetype) = 0;

    /**
     * AMIKLUSetOrdering sets the ordering for the sparse linear solver of the
     * forward problem
     *
     * @param[in] ordering ordering algorithm to reduce fill 0:AMD 1:COLAMD 2:
     * natural ordering @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMIKLUSetOrdering(int ordering) = 0;

    /**
     * AMIKLUSetOrderingB sets the ordering for the sparse linear solver of the
     * backward problem
     *
     * @param[in] which identifier of the backwards problem @type int
     * @param[in] ordering ordering algorithm to reduce fill 0:AMD 1:COLAMD 2:
     * natural ordering @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMIKLUSetOrderingB(int which, int ordering) = 0;

    /**
     * AMIKLUB attaches a sparse linear solver to the forward problem
     *
     * @param[in] which identifier of the backwards problem @type int
     * @param[in] nx number of state variables @type int
     * @param[in] nnz number of nonzero entries in the jacobian @type int
     * @param[in] sparsetype sparse storage type, CSC_MAT for column matrix,
     * CSR_MAT for row matrix @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMIKLUB(int which, int nx, int nnz, int sparsetype) = 0;

    /**
     * AMIGetNumSteps reports the number of solver steps
     *
     * @param[in] ami_mem pointer to the solver memory object (can be from
     * forward or backward problem) @type void
     * @param[out] numsteps output array @type long int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMIGetNumSteps(void *ami_mem, long int *numsteps) = 0;

    /**
     * AMIGetNumRhsEvals reports the number of right hand evaluations
     *
     * @param[in] ami_mem pointer to the solver memory object (can be from
     * forward or backward problem) @type void
     * @param[out] numrhsevals output array @type long int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMIGetNumRhsEvals(void *ami_mem, long int *numrhsevals) = 0;

    /**
     * AMIGetNumErrTestFails reports the number of local error test failures
     *
     * @param[in] ami_mem pointer to the solver memory object (can be from
     * forward or backward problem) @type void
     * @param[out] numerrtestfails output array @type long int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMIGetNumErrTestFails(void *ami_mem,
                                      long int *numerrtestfails) = 0;

    /**
     * AMIGetNumNonlinSolvConvFails reports the number of nonlinear convergence
     * failures
     *
     * @param[in] ami_mem pointer to the solver memory object (can be from
     * forward or backward problem) @type void
     * @param[out] numnonlinsolvconvfails output array @type long int
     * @return status flag indicating success of execution @type int
     */
    virtual int
    AMIGetNumNonlinSolvConvFails(void *ami_mem,
                                 long int *numnonlinsolvconvfails) = 0;

    /**
     * AMIGetLastOrder reports the order of the integration method during the
     * last internal step
     *
     * @param[in] ami_mem pointer to the solver memory object (can be from
     * forward or backward problem) @type void
     * @param[out] order output array @type long int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMIGetLastOrder(void *ami_mem, int *order) = 0;

    /**
     * AMIGetAdjBmem retrieves the solver memory object for the backward problem
     *
     * @param[in] which identifier of the backwards problem @type int
     * @param[in] ami_mem pointer to the forward solver memory object @type void
     * @return ami_memB pointer to the backward solver memory object @type void
     */
    virtual void *AMIGetAdjBmem(void *ami_mem, int which) = 0;

    int setLinearSolver(const UserData *udata, Model *model);

    /** pointer to ami memory block */
    void *ami_mem = nullptr;
};

#endif // AMICISOLVER_H
