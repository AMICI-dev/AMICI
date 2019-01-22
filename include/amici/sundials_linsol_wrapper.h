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
    SUNLinSolWrapper& operator=(const SUNLinSolWrapper& other);

    /**
     * @brief Move assignment
     * @param other
     * @return
     */
    SUNLinSolWrapper& operator=(SUNLinSolWrapper&& other) noexcept;

    /**
     * @brief Get the wrapped SUNLinSol
     * @return SlsMat
     */
    SUNLinearSolver get() const;

    /**
     * @brief getType
     * @return
     */
    SUNLinearSolver_Type getType() const;

    /**
     * @brief initialize
     * @return
     */
    int initialize();

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

} // namespace amici
#endif // AMICI_SUNDIALS_LINSOL_WRAPPER_H

