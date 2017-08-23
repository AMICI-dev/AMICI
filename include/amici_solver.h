#ifndef AMICISOLVER_H
#define AMICISOLVER_H

#include <nvector/nvector_serial.h> // DlsMat
#include <sundials/sundials_sparse.h> // SlsMat

class ReturnData;
class UserData;
class TempData;
class Model;

class Solver
{
public:
    Solver() {

    }

    virtual ~Solver();

    /**
     * @brief setupAMIs initialises the ami memory object
     * @param[out] status flag indicating success of execution @type int
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[in] tdata pointer to the temporary data struct @type TempData
     * @return ami_mem pointer to the cvodes/idas memory block
     */

    int setupAMI(UserData *udata, TempData *tdata, Model *model);

    /**
     * setupAMIB initialises the AMI memory object for the backwards problem
     * @param[out] status flag indicating success of execution @type int
     * @param[in] ami_mem pointer to the solver memory object of the forward problem
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[in] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */

    int setupAMIB(UserData *udata, TempData *tdata, Model *model);

    virtual int AMIGetSens(realtype *tret, N_Vector *yySout) = 0;

    /**
     * getDiagnosis extracts diagnosis information from solver memory block and writes them into the return data struct
     *
     * @param[out]
     * @param[in] it time-point index @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @return status flag indicating success of execution @type int
     */

    int getDiagnosis(int it, ReturnData *rdata);

    /**
     * getDiagnosisB extracts diagnosis information from solver memory block and writes them into the return data struct for the backward problem
     *
     * @param[in] it time-point index @type int
     * @param[in] ami_mem pointer to the solver memory block @type *void
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[out] rdata pointer to the return data struct @type ReturnData
     * @param[out] tdata pointer to the temporary data struct @type TempData
     * @return status flag indicating success of execution @type int
     */

    int getDiagnosisB(int it, ReturnData *rdata, TempData *tdata);

    virtual int AMIGetRootInfo(int *rootsfound) = 0;

    virtual int AMIReInit(realtype t0, N_Vector yy0, N_Vector yp0) = 0;

    virtual int AMISensReInit(int ism, N_Vector *yS0, N_Vector *ypS0) = 0;

    virtual int AMICalcIC(realtype tout1) = 0;

    virtual int AMICalcICB(int which, realtype tout1, N_Vector xB, N_Vector dxB) = 0;

    virtual int AMISolve(realtype tout, N_Vector yret, N_Vector ypret, realtype *tret, int itask) = 0;

    virtual int AMISolveF(realtype tout, N_Vector yret, N_Vector ypret, realtype *tret, int itask, int *ncheckPtr) = 0;

    virtual int AMISolveB(realtype tBout, int itaskB) = 0;

    virtual int AMISetStopTime(realtype tstop) = 0;

//    virtual int AMIRootInit(int nrtfn, RootFn ptr) = 0;

    virtual int AMIReInitB(int which, realtype tB0, N_Vector yyB0, N_Vector ypB0) = 0;

    virtual int AMIGetB(int which, realtype *tret, N_Vector yy, N_Vector yp) = 0;

    virtual int AMIGetQuadB(int which, realtype *tret, N_Vector qB) = 0;

    virtual int AMIQuadReInitB(int which, N_Vector yQB0) = 0;

    virtual int turnOffRootFinding() = 0;

protected:
    virtual int wrap_init(N_Vector x, N_Vector dx, realtype t) = 0;

    // TODO: check if model has adjoint sensitivities, else return -1
    virtual int wrap_binit(int which, N_Vector xB, N_Vector dxB, realtype t) = 0;

    // TODO: check if model has adjoint sensitivities, else return -1
    virtual int wrap_qbinit(int which, N_Vector qBdot) = 0;

    virtual int wrap_RootInit(int ne) = 0;

    // TODO: check if model has forward sensitivities, else return -1
    virtual int wrap_SensInit1(N_Vector *sx, N_Vector *sdx, UserData *udata) = 0;

    virtual int wrap_SetDenseJacFn() = 0;

    virtual int wrap_SetSparseJacFn() = 0;

    virtual int wrap_SetBandJacFn() = 0;

    virtual int wrap_SetJacTimesVecFn() = 0;

    // TODO: check if model has adjoint sensitivities, else return -1
    virtual int wrap_SetDenseJacFnB(int which) = 0;

    // TODO: check if model has adjoint sensitivities, else return -1
    virtual int wrap_SetSparseJacFnB(int which) = 0;

    // TODO: check if model has adjoint sensitivities, else return -1
    virtual int wrap_SetBandJacFnB(int which) = 0;

    // TODO: check if model has adjoint sensitivities, else return -1
    virtual int wrap_SetJacTimesVecFnB(int which) = 0;

    static void wrap_ErrHandlerFn(int error_code, const char *module, const char *function, char *msg, void *eh_data);

    virtual void *AMICreate(int lmm, int iter) = 0;

    virtual int AMISStolerances(double rtol,double atol) = 0;

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

    virtual int AMICreateB(int lmm, int iter, int *which) = 0;

    virtual int AMISStolerancesB(int which, realtype relTolB, realtype absTolB) = 0;

    virtual int AMIQuadSStolerancesB(int which, realtype reltolQB, realtype abstolQB) = 0;

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

    virtual int AMIGetNumErrTestFails(void *ami_mem, long int *numerrtestfails) = 0;

    virtual int AMIGetNumNonlinSolvConvFails(void *ami_mem, long int *numnonlinsolvconvfails) = 0;

    virtual int AMIGetLastOrder(void *ami_mem, int *order) = 0;

    virtual void *AMIGetAdjBmem(void *ami_mem, int which) = 0;

    int setLinearSolver(const UserData *udata, Model *model);

    /** pointer to ami memory block */
    void *ami_mem = nullptr;

};



#endif // AMICISOLVER_H
