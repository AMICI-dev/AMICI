#ifndef AMICISOLVER_H
#define AMICISOLVER_H

#include <nvector/nvector_serial.h>
class ReturnData;
class UserData;
class TempData;
class Model;

// TODO: get out of here
typedef int (*RootFn)(realtype t, N_Vector y, realtype *gout, void *user_data);

#include <sundials/sundials_spgmr.h>
// TODO: don't use cvodes includes here
#include <cvodes/cvodes_spils.h>


class Solver
{
public:
    Solver() {

    }

    virtual int wrap_init(void *mem, N_Vector x, N_Vector dx, realtype t) = 0;

    // TODO: check if model has adjoint sensitivities, else return -1
    virtual int wrap_binit(void *mem, int which, N_Vector xB, N_Vector dxB, realtype t) = 0;

    // TODO: check if model has adjoint sensitivities, else return -1

    virtual int wrap_qbinit(void *mem, int which, N_Vector qBdot) = 0;

    virtual int wrap_RootInit(void *mem, UserData *udata) = 0;

    virtual int wrap_SensInit1(void *mem, N_Vector *sx, N_Vector *sdx, UserData *udata) = 0;

    /**
     * @brief setupAMIs initialises the ami memory object
     * @param[out] status flag indicating success of execution @type int
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[in] tdata pointer to the temporary data struct @type TempData
     * @return ami_mem pointer to the cvodes/idas memory block
     */

    void *setupAMI(UserData *udata, TempData *tdata);

    /**
     * setupAMIB initialises the AMI memory object for the backwards problem
     * @param[out] status flag indicating success of execution @type int
     * @param[in] ami_mem pointer to the solver memory object of the forward problem
     * @param[in] udata pointer to the user data struct @type UserData
     * @param[in] tdata pointer to the temporary data struct @type TempData
     * @return ami_mem pointer to the cvodes/idas memory block for the backward problem
     */

    int setupAMIB(void *ami_mem, UserData *udata, TempData *tdata);

    static void wrap_ErrHandlerFn(int error_code, const char *module, const char *function, char *msg, void *eh_data);

    virtual void *AMICreate(int lmm, int iter) = 0;

    virtual int AMISStolerances(void *mem, double rtol,double atol) = 0;

    virtual int AMISensEEtolerances(void *mem) = 0;

    virtual int AMISetSensErrCon(void *mem,bool error_corr) = 0;

    virtual int AMISetQuadErrConB(void *mem,int which, bool flag) = 0;

    virtual int AMIGetRootInfo(void *mem,int *rootsfound) = 0;

    virtual int AMISetErrHandlerFn(void *mem) = 0;

    virtual int AMISetUserData(void *mem, void *user_data) = 0;

    virtual int AMISetUserDataB(void *mem, int which, void *user_data) = 0;

    virtual int AMISetMaxNumSteps(void *mem, long int mxsteps) = 0;

    virtual int AMISetStabLimDet(void *mem, int stldet) = 0;

    virtual int AMISetStabLimDetB(void *mem, int which, int stldet) = 0;

    virtual int AMISetId(void *mem, N_Vector id) = 0;

    virtual int AMISetSuppressAlg(void *mem, bool flag) = 0;

    virtual int AMIReInit(void *mem, realtype t0, N_Vector yy0, N_Vector yp0) = 0;

    virtual int AMISensReInit(void *mem, int ism, N_Vector *yS0, N_Vector *ypS0) = 0;

    virtual int AMISetSensParams(void *mem, realtype *p, realtype *pbar, int *plist) = 0;

    virtual int AMIGetDky(void *mem, realtype t, int k, N_Vector dky) = 0;

    virtual int AMIGetSens(void *mem, realtype *tret, N_Vector *yySout) = 0;

    virtual int AMIRootInit(void *mem, int nrtfn, RootFn ptr) = 0;

    virtual void AMIFree(void **mem) = 0;

    virtual int AMIAdjInit(void *mem, long int steps, int interp) = 0;

    virtual int AMICreateB(void *mem, int lmm, int iter, int *which) = 0;

    virtual int AMIReInitB(void *mem, int which, realtype tB0, N_Vector yyB0, N_Vector ypB0) = 0;

    virtual int AMISStolerancesB(void *mem, int which, realtype relTolB, realtype absTolB) = 0;

    virtual int AMIQuadReInitB(void *mem, int which, N_Vector yQB0) = 0;

    virtual int AMIQuadSStolerancesB(void *mem, int which, realtype reltolQB, realtype abstolQB) = 0;

    virtual int AMISolve(void *mem, realtype tout, N_Vector yret, N_Vector ypret, realtype *tret, int itask) = 0;

    virtual int AMISolveF(void *mem, realtype tout, N_Vector yret, N_Vector ypret, realtype *tret, int itask, int *ncheckPtr) = 0;

    virtual int AMISolveB(void *mem, realtype tBout, int itaskB) = 0;

    virtual int AMISetMaxNumStepsB(void *mem, int which, long int mxstepsB) = 0;

    virtual int AMIGetB(void *mem, int which, realtype *tret, N_Vector yy, N_Vector yp) = 0;

    virtual int AMIGetQuadB(void *mem, int which, realtype *tret, N_Vector qB) = 0;

    virtual int AMIDense(void *mem, int nx) = 0;

    virtual int AMIDenseB(void *mem, int which, int nx) = 0;

    virtual int AMIBand(void *mem, int nx, int ubw, int lbw) = 0;

    virtual int AMIBandB(void *mem, int which, int nx, int ubw, int lbw) = 0;

    virtual int AMIDiag(void *mem) = 0;

    virtual int AMIDiagB(void *mem, int which) = 0;

    virtual int AMISpgmr(void *mem, int prectype, int maxl) = 0;

    virtual int AMISpgmrB(void *mem, int which, int prectype, int maxl) = 0;

    virtual int AMISpbcg(void *mem, int prectype, int maxl) = 0;

    virtual int AMISpbcgB(void *mem, int which, int prectype, int maxl) = 0;

    virtual int AMISptfqmr(void *mem, int prectype, int maxl) = 0;

    virtual int AMISptfqmrB(void *mem, int which, int prectype, int maxl) = 0;

    virtual int AMIKLU(void *mem, int nx, int nnz, int sparsetype) = 0;

    virtual int AMIKLUSetOrdering(void *mem, int ordering) = 0;

    virtual int AMIKLUSetOrderingB(void *mem, int which, int ordering) = 0;

    virtual int AMIKLUB(void *mem, int which, int nx, int nnz, int sparsetype) = 0;

    virtual int AMIGetNumSteps(void *mem, long int *numsteps) = 0;

    virtual int AMIGetNumRhsEvals(void *mem, long int *numrhsevals) = 0;

    virtual int AMIGetNumErrTestFails(void *mem, long int *numerrtestfails) = 0;

    virtual int AMIGetNumNonlinSolvConvFails(void *mem, long int *numnonlinsolvconvfails) = 0;

    virtual int AMIGetLastOrder(void *mem,int *order) = 0;

    virtual void *AMIGetAdjBmem(void *mem, int which) = 0;

    virtual int AMICalcIC(void *mem, realtype tout1) = 0;

    virtual int AMICalcICB(void *mem, int which, realtype tout1, N_Vector xB, N_Vector dxB) = 0;

    virtual int AMISetStopTime(void *mem, realtype tstop) = 0;

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

    int getDiagnosis(int it, void *ami_mem, ReturnData *rdata);

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

    int getDiagnosisB(int it, void *ami_mem, UserData *udata, ReturnData *rdata, TempData *tdata);

    virtual ~Solver() {

    }
};



#endif // AMICISOLVER_H
