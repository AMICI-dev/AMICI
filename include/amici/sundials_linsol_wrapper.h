#ifndef AMICI_SUNDIALS_LINSOL_WRAPPER_H
#define AMICI_SUNDIALS_LINSOL_WRAPPER_H

#include "amici/exception.h"
#include "amici/sundials_matrix_wrapper.h"

#include <sundials/sundials_linearsolver.h> // SUNLinearSolver
#include <sunlinsol/sunlinsol_band.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_klu.h>
#include <sunlinsol/sunlinsol_pcg.h>
#include <sunlinsol/sunlinsol_spbcgs.h>
#include <sunlinsol/sunlinsol_spfgmr.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunlinsol/sunlinsol_sptfqmr.h>

#include <sundials/sundials_nonlinearsolver.h>
#include<sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include<sunnonlinsol/sunnonlinsol_newton.h>


namespace amici {

/**
 * @brief A RAII wrapper for SUNLinearSolver structs.
 *
 * More members to be added as required.
 * Start accessing through SUNLinSolWrapper::linsol().
 *
 * TODO: check return values
 */
class SUNLinSolWrapper {
public:
    /**
     * @brief SUNLinSolWrapper from existing SUNLinearSolver
     * @param linsol_
     */
    explicit SUNLinSolWrapper(SUNLinearSolver linsol_);

    virtual ~SUNLinSolWrapper();

    /**
     * @brief Copy constructor
     * @param other
     */
    SUNLinSolWrapper(const SUNLinSolWrapper& other) = delete;

    /**
     * @brief Move constructor
     * @param other
     */
    SUNLinSolWrapper(SUNLinSolWrapper&& other) noexcept;

    /**
     * @brief Copy assignment
     * @param other
     * @return
     */
    SUNLinSolWrapper& operator=(const SUNLinSolWrapper& other) = delete;

    /**
     * @brief Move assignment
     * @param other
     * @return
     */
    SUNLinSolWrapper& operator=(SUNLinSolWrapper&& other) noexcept;

    /**
     * @brief Get the wrapped SUNLinSol
     * @return SUNLinearSolver
     */
    SUNLinearSolver get() const;

    /**
     * @brief getType
     * @return
     */
    SUNLinearSolver_Type getType() const;

    /**
     * @brief setup
     * @param A
     * @return
     */
    int setup(SUNMatrix A);

    /**
     * @brief setup
     * @param A
     * @return
     */
    int setup(SUNMatrixWrapper A);

    /**
     * @brief Solve
     * @param A
     * @param x
     * @param b
     * @param tol
     * @return
     */
    int Solve(SUNMatrix A, N_Vector x, N_Vector b, realtype tol);

    /**
     * @brief getLastFlag
     * @return
     */
    long int getLastFlag();

    /**
     * @brief space_Dense
     * @param lenrwLS
     * @param leniwLS
     * @return
     */
    int space(long int *lenrwLS, long int *leniwLS);

protected:

    /**
     * @brief initialize
     * @return
     */
    int initialize();

    /** the wrapper solver */
    SUNLinearSolver linsol = nullptr;
};


class SUNLinSolBand: public SUNLinSolWrapper {
    /**
     * @brief SUNLinSolBand
     * @param y
     * @param A
     */
    SUNLinSolBand(N_Vector y, SUNMatrix A)
        : SUNLinSolWrapper (SUNLinSol_Band(y, A))
    {
        if(!linsol)
            throw AmiException("Failed to create solver.");
        initialize();
    }
};


class SUNLinSolDense: public SUNLinSolWrapper {
    /**
     * @brief SUNLinSolDense
     * @param y
     * @param A
     */
    SUNLinSolDense(N_Vector y, SUNMatrix A)
        : SUNLinSolWrapper (SUNLinSol_Dense(y, A))
    {
        if(!linsol)
            throw AmiException("Failed to create solver.");
        initialize();
    }
};


class SUNLinSolKLU : public SUNLinSolWrapper {
    /**
     * @brief SUNLinSolKLU
     * @param y
     * @param A
     */
    SUNLinSolKLU(N_Vector y, SUNMatrix A)
        : SUNLinSolWrapper(SUNLinSol_KLU(y, A))
    {
        if(!linsol)
            throw AmiException("Failed to create solver.");
        initialize();
    }
};


class SUNLinSolPCG: public SUNLinSolWrapper {
    /**
     * @brief SUNLinSolPCG
     * @param y
     * @param pretype
     * @param maxl
     */
    SUNLinSolPCG(N_Vector y, int pretype, int maxl)
        : SUNLinSolWrapper(SUNLinSol_PCG(y, pretype, maxl))
    {
        if(!linsol)
            throw AmiException("Failed to create solver.");
        initialize();
    }

    int setATimes(void* A_data, ATimesFn ATimes)
    {
        return SUNLinSolSetATimes_PCG(linsol, A_data, ATimes);
    }

    int setPreconditioner(void* P_data, PSetupFn Pset, PSolveFn Psol) {
        return SUNLinSolSetPreconditioner_PCG(linsol, P_data, Pset, Psol);
    }

    int setScalingVectors(N_Vector s, N_Vector nul) {
        return SUNLinSolSetScalingVectors_PCG(linsol, s, nul);
    }

    int getNumIters() {
        return SUNLinSolNumIters_PCG(linsol);
    }

    realtype getResNorm() {
        return SUNLinSolResNorm_PCG(linsol);
    }

    N_Vector getResid() {
        return SUNLinSolResid_PCG(linsol);
    }

};


class SUNLinSolSPBCGS : public SUNLinSolWrapper {
    SUNLinSolSPBCGS(N_Vector y, int pretype, int maxl)
        : SUNLinSolWrapper(SUNLinSol_SPBCGS(y, pretype, maxl))
    {
        if(!linsol)
            throw AmiException("Failed to create solver.");
        initialize();
    }

    int setATimes(void* A_data, ATimesFn ATimes)
    {
        return SUNLinSolSetATimes_SPBCGS(linsol, A_data, ATimes);
    }

    int setPreconditioner(void* P_data, PSetupFn Pset, PSolveFn Psol) {
        return SUNLinSolSetPreconditioner_SPBCGS(linsol, P_data, Pset, Psol);
    }

    int setScalingVectors(N_Vector s, N_Vector nul) {
        return SUNLinSolSetScalingVectors_SPBCGS(linsol, s, nul);
    }

    int getNumIters() {
        return SUNLinSolNumIters_SPBCGS(linsol);
    }

    realtype getResNorm() {
        return SUNLinSolResNorm_SPBCGS(linsol);
    }

    N_Vector getResid() {
        return SUNLinSolResid_SPBCGS(linsol);
    }

};


class SUNLinSolSPFGMR: public SUNLinSolWrapper {
    /**
     * @brief SUNLinSolSPFGMR
     * @param y
     * @param pretype
     * @param maxl
     */
    SUNLinSolSPFGMR(N_Vector y, int pretype, int maxl)
        : SUNLinSolWrapper(SUNLinSol_SPFGMR(y, pretype, maxl))
    {
        if(!linsol)
            throw AmiException("Failed to create solver.");
        initialize();
    }

    int setATimes(void* A_data, ATimesFn ATimes)
    {
        return SUNLinSolSetATimes_SPFGMR(linsol, A_data, ATimes);
    }

    int setPreconditioner(void* P_data, PSetupFn Pset, PSolveFn Psol) {
        return SUNLinSolSetPreconditioner_SPFGMR(linsol, P_data, Pset, Psol);
    }

    int setScalingVectors(N_Vector s, N_Vector nul) {
        return SUNLinSolSetScalingVectors_SPFGMR(linsol, s, nul);
    }

    int getNumIters() {
        return SUNLinSolNumIters_SPFGMR(linsol);
    }

    realtype getResNorm() {
        return SUNLinSolResNorm_SPFGMR(linsol);
    }

    N_Vector getResid() {
        return SUNLinSolResid_SPFGMR(linsol);
    }

};


class SUNLinSolSPGMR: public SUNLinSolWrapper {
    /**
     * @brief SUNLinSolSPGMR
     * @param y
     * @param pretype
     * @param maxl
     */
    SUNLinSolSPGMR(N_Vector y, int pretype, int maxl)
        : SUNLinSolWrapper(SUNLinSol_SPGMR(y, pretype, maxl))
    {
        if(!linsol)
            throw AmiException("Failed to create solver.");
        initialize();
    }

    int setATimes(void* A_data, ATimesFn ATimes)
    {
        return SUNLinSolSetATimes_SPGMR(linsol, A_data, ATimes);
    }

    int setPreconditioner(void* P_data, PSetupFn Pset, PSolveFn Psol) {
        return SUNLinSolSetPreconditioner_SPGMR(linsol, P_data, Pset, Psol);
    }

    int setScalingVectors(N_Vector s, N_Vector nul) {
        return SUNLinSolSetScalingVectors_SPGMR(linsol, s, nul);
    }

    int getNumIters() {
        return SUNLinSolNumIters_SPGMR(linsol);
    }

    realtype getResNorm() {
        return SUNLinSolResNorm_SPGMR(linsol);
    }

    N_Vector getResid() {
        return SUNLinSolResid_SPGMR(linsol);
    }
};


class SUNLinSolSPTFQMR: public SUNLinSolWrapper {
    /**
     * @brief SUNLinSolSPTFQMR
     * @param y
     * @param pretype
     * @param maxl
     */
    SUNLinSolSPTFQMR(N_Vector y, int pretype, int maxl)
        : SUNLinSolWrapper(SUNLinSol_SPTFQMR(y, pretype, maxl))
    {
        if(!linsol)
            throw AmiException("Failed to create solver.");
        initialize();
    }

    int setATimes(void* A_data, ATimesFn ATimes)
    {
        return SUNLinSolSetATimes_SPTFQMR(linsol, A_data, ATimes);
    }

    int setPreconditioner(void* P_data, PSetupFn Pset, PSolveFn Psol) {
        return SUNLinSolSetPreconditioner_SPTFQMR(linsol, P_data, Pset, Psol);
    }

    int setScalingVectors(N_Vector s, N_Vector nul) {
        return SUNLinSolSetScalingVectors_SPTFQMR(linsol, s, nul);
    }

    int getNumIters() {
        return SUNLinSolNumIters_SPTFQMR(linsol);
    }

    realtype getResNorm() {
        return SUNLinSolResNorm_SPTFQMR(linsol);
    }

    N_Vector getResid() {
        return SUNLinSolResid_SPTFQMR(linsol);
    }

};



/**
 * @brief A RAII wrapper for SUNNonLinearSolver structs.
 *
 * TODO: check return values
 */
class SUNNonLinSolWrapper {
public:
    /**
     * @brief SUNNonLinSolWrapper from existing SUNNonlinearSolver
     * @param linsol_
     */
    explicit SUNNonLinSolWrapper(SUNNonlinearSolver sol_);

    virtual ~SUNNonLinSolWrapper();

    /**
     * @brief Copy constructor
     * @param other
     */
    SUNNonLinSolWrapper(const SUNNonLinSolWrapper& other) = delete;

    /**
     * @brief Move constructor
     * @param other
     */
    SUNNonLinSolWrapper(SUNNonLinSolWrapper&& other) noexcept;

    /**
     * @brief Copy assignment
     * @param other
     * @return
     */
    SUNNonLinSolWrapper& operator=(const SUNNonLinSolWrapper& other) = delete;

    /**
     * @brief Move assignment
     * @param other
     * @return
     */
    SUNNonLinSolWrapper& operator=(SUNNonLinSolWrapper&& other) noexcept;

    /**
     * @brief Get the wrapped SUNNonlinearSolver
     * @return SUNNonlinearSolver
     */
    SUNNonlinearSolver get() const;

    /**
     * @brief getType
     * @return
     */
    SUNNonlinearSolver_Type getType() const;

    /**
     * @brief setup
     * @param y
     * @param mem
     * @return
     */
    int setup(N_Vector y, void* mem);


    int Solve(N_Vector y0, N_Vector y,
              N_Vector w, realtype tol,
              booleantype callLSetup, void *mem);

    int setSysFn(SUNNonlinSolSysFn SysFn) {
        return SUNNonlinSolSetSysFn(solver, SysFn);
    }

    int setLSetupFn(SUNNonlinSolLSetupFn SetupFn) {
        return SUNNonlinSolSetLSetupFn(solver, SetupFn);
    }

    int setLSolveFn(SUNNonlinSolLSolveFn SolveFn) {
        return SUNNonlinSolSetLSolveFn(solver, SolveFn);
    }

    int setConvTestFn(SUNNonlinSolConvTestFn CTestFn) {
        return SUNNonlinSolSetConvTestFn(solver, CTestFn);
    }

    int setMaxIters(int maxiters) {
        return SUNNonlinSolSetMaxIters(solver, maxiters);
    }

    int getNumIters(long int *niters) {
        return SUNNonlinSolGetNumIters(solver, niters);
    }

    int getCurIter(int *iter) {
        return SUNNonlinSolGetCurIter(solver, iter);
    }

    int getNumConvFails(long int *nconvfails) {
        return SUNNonlinSolGetNumConvFails(solver, nconvfails);
    }

protected:

    /**
     * @brief initialize
     * @return
     */
    int initialize();

    /** the wrapper solver */
    SUNNonlinearSolver solver = nullptr;
};



class SUNNonLinSolNewton: public SUNNonLinSolWrapper {
    SUNNonLinSolNewton(N_Vector y)
        :SUNNonLinSolWrapper(SUNNonlinSol_Newton(y))
    {
        initialize();
    }

    SUNNonLinSolNewton(int count, N_Vector y)
        :SUNNonLinSolWrapper(SUNNonlinSol_NewtonSens(count, y))
    {
        initialize();
    }

    int getSysFn(SUNNonlinSolSysFn *SysFn) {
        return SUNNonlinSolGetSysFn_Newton(solver, SysFn);
    }
};


class SUNNonLinSolFixedPoint: public SUNNonLinSolWrapper {

    SUNNonLinSolFixedPoint(N_Vector y, int m)
        :SUNNonLinSolWrapper(SUNNonlinSol_FixedPoint(y, m))
    {
        initialize();
    }

    SUNNonLinSolFixedPoint(int count, N_Vector y, int m)
        :SUNNonLinSolWrapper(SUNNonlinSol_FixedPointSens(count, y, m))
    {
        initialize();
    }

    int getSysFn(SUNNonlinSolSysFn *SysFn) {
        return SUNNonlinSolGetSysFn_FixedPoint(solver, SysFn);
    }
};


} // namespace amici
#endif // AMICI_SUNDIALS_LINSOL_WRAPPER_H

