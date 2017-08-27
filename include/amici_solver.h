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
  provides a generic interface to CVode and IDA solvers, individual realizations are realized in the CVodeSolver and the IDASolver class.
*/
class Solver
{
public:
    Solver() {}

    virtual ~Solver();

    int setupAMI(UserData *udata, TempData *tdata, Model *model);

    int setupAMIB(UserData *udata, TempData *tdata, Model *model);

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
     * @param[out] rootsfound array with flags indicating whether the respective event occured
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
     * AMISensReInit reinitializes the state sensitivites in the solver after an event occurence
     *
     * @param[in] ism sensitivity mode @type realtype
     * @param[in] yS0 new state sensitivity @type N_Vector
     * @param[in] ypS0 new derivative state sensitivities (DAE only) @type N_Vector
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISensReInit(int ism, N_Vector *yS0, N_Vector *ypS0) = 0;

   /**
     * AMICalcIC calculates consistent initial conditions, assumes initial states to be correct (DAE only)
     *
     * @param[in] tout1 next timepoint to be computed (sets timescale) @type realtype
     * @return status flag indicating success of execution @type int
     */
    virtual int AMICalcIC(realtype tout1) = 0;

   /**
     * AMICalcIBC calculates consistent initial conditions for the backwards problem, assumes initial states to be correct (DAE only)
     *
     * @param[in] which identifier of the backwards problem @type int
     * @param[in] tout1 next timepoint to be computed (sets timescale) @type realtype
     * @param[in] xB states of final solution of the forward problem @type N_Vector
     * @param[in] dxB derivative states of final solution of the forward problem (DAE only) @type N_Vector
     * @return status flag indicating success of execution @type int
     */
    virtual int AMICalcICB(int which, realtype tout1, N_Vector xB,
                           N_Vector dxB) = 0;

   /**
     * AMISolve solves the forward problem until a predefined timepoint
     *
     * @param[in] tout timepoint until which simulation should be performed @type realtype
     * @param[in] yret states @type N_Vector
     * @param[in] ypret derivative states (DAE only) @type N_Vector
     * @param[in,out] tret pointer to the time variable @type realtype
     * @param[in] itask task identifier, can be CV_NORMAL or CV_ONE_STEP @type realtype
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISolve(realtype tout, N_Vector yret, N_Vector ypret,
                         realtype *tret, int itask) = 0;

   /**
     * AMISolveF solves the forward problem until a predefined timepoint (adjoint only)
     *
     * @param[in] tout timepoint until which simulation should be performed @type realtype
     * @param[in] yret states @type N_Vector
     * @param[in] ypret derivative states (DAE only) @type N_Vector
     * @param[in,out] tret pointer to the time variable @type realtype
     * @param[in] itask task identifier, can be CV_NORMAL or CV_ONE_STEP @type realtype
     * @param[in] ncheckPtr pointer to a number that counts the internal checkpoints @type realtype
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISolveF(realtype tout, N_Vector yret, N_Vector ypret,
                          realtype *tret, int itask, int *ncheckPtr) = 0;

   /**
     * AMISolveB solves the backward problem until a predefined timepoint (adjoint only)
     *
     * @param[in] tBout timepoint until which simulation should be performed @type realtype
     * @param[in] itaskB task identifier, can be CV_NORMAL or CV_ONE_STEP @type realtype
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISolveB(realtype tBout, int itaskB) = 0;

   /**
     * AMISetStopTime sets a timepoint at which the simulation will be stopped
     *
     * @param[in] tstop timepoint until which simulation should be performed @type realtype
     * @return status flag indicating success of execution @type int
     */
    virtual int AMISetStopTime(realtype tstop) = 0;

    //    virtual int AMIRootInit(int nrtfn, RootFn ptr) = 0;

   /**
     * AMIReInitB reinitializes the adjoint states after an event occurence
     *
     * @param[in] which identifier of the backwards problem @type int
     * @param[in] t0 new timepoint @type realtype
     * @param[in] yyB0 new adjoint state variables @type N_Vector
     * @param[in] ypB0 new adjoint derivative state variables (DAE only) @type N_Vector
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
     * @param[in] yp adjoint derivative state variables (DAE only) @type N_Vector
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
     * wrap_init initialises the states at the specified initial timepoint
     *
     * @param[in] x initial state variables @type N_Vector
     * @param[in] dx initial derivative state variables (DAE only) @type N_Vector
     * @param[in] t initial timepoint @type realtype
     * @return status flag indicating success of execution @type int
     */
    virtual int wrap_init(N_Vector x, N_Vector dx, realtype t) = 0;

    // TODO: check if model has adjoint sensitivities, else return -1

    /**
     * wrap_binit initialises the adjoint states at the specified final timepoint
     *
     * @param[in] which identifier of the backwards problem @type int
     * @param[in] xB initial adjoint state variables @type N_Vector
     * @param[in] dxB initial adjoint derivative state variables (DAE only) @type N_Vector
     * @param[in] t final timepoint @type realtype
     * @return status flag indicating success of execution @type int
     */
    virtual int wrap_binit(int which, N_Vector xB, N_Vector dxB,
                           realtype t) = 0;

    // TODO: check if model has adjoint sensitivities, else return -1

    /**
     * wrap_qbinit initialises the quadrature states at the specified final timepoint
     *
     * @param[in] which identifier of the backwards problem @type int
     * @param[in] qBdot initial adjoint quadrature state variables @type N_Vector
     * @return status flag indicating success of execution @type int
     */
    virtual int wrap_qbinit(int which, N_Vector qBdot) = 0;

    /**
     * wrap_RootInit initialises the rootfinding for events
     *
     * @param[in] ne number of different events @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int wrap_RootInit(int ne) = 0;

    // TODO: check if model has forward sensitivities, else return -1
    /**
     * wrap_SensInit1 initialises the sensitivities at the specified initial timepoint
     *
     * @param[in] sx initial state sensitivities @type N_Vector
     * @param[in] dx initial derivative state sensitivities (DAE only) @type N_Vector
     * @param[in] udata initial derivative state sensitivities (DAE only) @type N_Vector
     * @return status flag indicating success of execution @type int
     */
    virtual int wrap_SensInit1(N_Vector *sx, N_Vector *sdx,
                               const UserData *udata) = 0;

    /**
     * wrap_SetDenseJacFn sets the dense Jacobian function
     *
     * @return status flag indicating success of execution @type int
     */
    virtual int wrap_SetDenseJacFn() = 0;

    /**
     * wrap_SetSparseJacFn sets the sparse Jacobian function
     *
     * @return status flag indicating success of execution @type int
     */
    virtual int wrap_SetSparseJacFn() = 0;

    /**
     * wrap_SetBandJacFn sets the banded Jacobian function
     *
     * @return status flag indicating success of execution @type int
     */
    virtual int wrap_SetBandJacFn() = 0;

    /**
     * wrap_SetJacTimesVecFn sets the Jacobian vector multiplication function
     *
     * @return status flag indicating success of execution @type int
     */
    virtual int wrap_SetJacTimesVecFn() = 0;

    // TODO: check if model has adjoint sensitivities, else return -1
    /**
     * wrap_SetDenseJacFn sets the dense Jacobian function
     *
     * @param[in] which identifier of the backwards problem @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int wrap_SetDenseJacFnB(int which) = 0;

    // TODO: check if model has adjoint sensitivities, else return -1
    /**
     * wrap_SetSparseJacFn sets the sparse Jacobian function
     *
     * @param[in] which identifier of the backwards problem @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int wrap_SetSparseJacFnB(int which) = 0;

    // TODO: check if model has adjoint sensitivities, else return -1
    /**
     * wrap_SetBandJacFn sets the banded Jacobian function
     *
     * @param[in] which identifier of the backwards problem @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int wrap_SetBandJacFnB(int which) = 0;

    // TODO: check if model has adjoint sensitivities, else return -1
    /**
     * wrap_SetJacTimesVecFn sets the Jacobian vector multiplication function
     *
     * @param[in] which identifier of the backwards problem @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int wrap_SetJacTimesVecFnB(int which) = 0;

    static void wrap_ErrHandlerFn(int error_code, const char *module,
                                  const char *function, char *msg,
                                  void *eh_data);

    /**
     * AMICreate specifies solver method and initializes solver memory
     *
     * @param[in] lmm linear multistep method CV_ADAMS or CV_BDF @type int
     * @param[in] lmm nonlinear solver method CV_NEWTON or CV_FUNCTIONAL @type int
     * @return status flag indicating success of execution @type int
     */
    virtual void *AMICreate(int lmm, int iter) = 0;

    virtual int AMISStolerances(double rtol, double atol) = 0;

    virtual int AMISensEEtolerances() = 0;

    virtual int AMISetSensErrCon(bool error_corr) = 0;

    virtual int AMISetQuadErrConB(int which, bool flag) = 0;

    virtual int AMISetErrHandlerFn() = 0;

    virtual int AMISetUserData(void *user_data) = 0;

    virtual int AMISetUserDataB(int which, void *user_data) = 0;

    virtual int AMISetMaxNumSteps(long int mxsteps) = 0;

    virtual int AMISetStabLimDet(int stldet) = 0;

    virtual int AMISetStabLimDetB(int which, int stldet) = 0;

    virtual int AMISetId(Model *model) = 0;

    virtual int AMISetSuppressAlg(bool flag) = 0;

    virtual int AMISetSensParams(realtype *p, realtype *pbar, int *plist) = 0;

    virtual int AMIGetDky(realtype t, int k, N_Vector dky) = 0;

    virtual void AMIFree() = 0;

    virtual int AMIAdjInit(long int steps, int interp) = 0;

    /**
     * AMICreateB specifies solver method and initializes solver memory
     *
     * @param[in] which identifier of the backwards problem @type int
     * @param[in] lmm linear multistep method CV_ADAMS or CV_BDF @type int
     * @param[in] lmm nonlinear solver method CV_NEWTON or CV_FUNCTIONAL @type int
     * @return status flag indicating success of execution @type int
     */
    virtual int AMICreateB(int lmm, int iter, int *which) = 0;

    virtual int AMISStolerancesB(int which, realtype relTolB,
                                 realtype absTolB) = 0;

    virtual int AMIQuadSStolerancesB(int which, realtype reltolQB,
                                     realtype abstolQB) = 0;

    virtual int AMISetMaxNumStepsB(int which, long int mxstepsB) = 0;

    virtual int AMIDense(int nx) = 0;

    virtual int AMIDenseB(int which, int nx) = 0;

    virtual int AMIBand(int nx, int ubw, int lbw) = 0;

    virtual int AMIBandB(int which, int nx, int ubw, int lbw) = 0;

    virtual int AMIDiag() = 0;

    virtual int AMIDiagB(int which) = 0;

    virtual int AMISpgmr(int prectype, int maxl) = 0;

    virtual int AMISpgmrB(int which, int prectype, int maxl) = 0;

    virtual int AMISpbcg(int prectype, int maxl) = 0;

    virtual int AMISpbcgB(int which, int prectype, int maxl) = 0;

    virtual int AMISptfqmr(int prectype, int maxl) = 0;

    virtual int AMISptfqmrB(int which, int prectype, int maxl) = 0;

    virtual int AMIKLU(int nx, int nnz, int sparsetype) = 0;

    virtual int AMIKLUSetOrdering(int ordering) = 0;

    virtual int AMIKLUSetOrderingB(int which, int ordering) = 0;

    virtual int AMIKLUB(int which, int nx, int nnz, int sparsetype) = 0;

    virtual int AMIGetNumSteps(void *ami_mem, long int *numsteps) = 0;

    virtual int AMIGetNumRhsEvals(void *ami_mem, long int *numrhsevals) = 0;

    virtual int AMIGetNumErrTestFails(void *ami_mem,
                                      long int *numerrtestfails) = 0;

    virtual int
    AMIGetNumNonlinSolvConvFails(void *ami_mem,
                                 long int *numnonlinsolvconvfails) = 0;

    virtual int AMIGetLastOrder(void *ami_mem, int *order) = 0;

    virtual void *AMIGetAdjBmem(void *ami_mem, int which) = 0;

    int setLinearSolver(const UserData *udata, Model *model);

    /** pointer to ami memory block */
    void *ami_mem = nullptr;
};

#endif // AMICISOLVER_H
