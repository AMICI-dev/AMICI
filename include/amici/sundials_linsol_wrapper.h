#ifndef AMICI_SUNDIALS_LINSOL_WRAPPER_H
#define AMICI_SUNDIALS_LINSOL_WRAPPER_H

#include "amici/exception.h"
#include "amici/sundials_matrix_wrapper.h"
#include "amici/vector.h"

#include <sunlinsol/sunlinsol_band.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_klu.h>
#include <sunlinsol/sunlinsol_pcg.h>
#include <sunlinsol/sunlinsol_spbcgs.h>
#include <sunlinsol/sunlinsol_spfgmr.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunlinsol/sunlinsol_sptfqmr.h>

#include<sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include<sunnonlinsol/sunnonlinsol_newton.h>


namespace amici {

/**
 * @brief A RAII wrapper for SUNLinearSolver structs.
 *
 * TODO: check return values
 */
class SUNLinSolWrapper {
public:
    SUNLinSolWrapper() = default;

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

    virtual SUNMatrix getMatrix() const;

protected:

    /**
     * @brief initialize
     * @return
     */
    int initialize();

    /** the wrapped solver */
    SUNLinearSolver linsol = nullptr;
};


class SUNLinSolBand: public SUNLinSolWrapper {
public:
    /**
     * @brief SUNLinSolBand
     * @param y
     * @param A
     */
    SUNLinSolBand(N_Vector y, SUNMatrix A);

    SUNLinSolBand(AmiVector const& x, int ubw, int lbw);

    SUNMatrix getMatrix() const override;

private:
    SUNMatrixWrapper A;
};


class SUNLinSolDense: public SUNLinSolWrapper {
public:
    SUNLinSolDense(AmiVector const& x);

    SUNMatrix getMatrix() const override;

private:
    SUNMatrixWrapper A;
};


class SUNLinSolKLU : public SUNLinSolWrapper {
public:
    /**
     * @brief SUNLinSolKLU
     * @param y
     * @param A
     */
    SUNLinSolKLU(N_Vector y, SUNMatrix A);

    SUNLinSolKLU(AmiVector const& x, int nnz, int sparsetype, StateOrdering ordering);

    SUNMatrix getMatrix() const override;

    void reInit(int nnz, int reinit_type);

    void setOrdering(StateOrdering ordering);

private:
    SUNMatrixWrapper A;
};


class SUNLinSolPCG: public SUNLinSolWrapper {
public:
    /**
     * @brief SUNLinSolPCG
     * @param y
     * @param pretype
     * @param maxl
     */
    SUNLinSolPCG(N_Vector y, int pretype, int maxl);

    int setATimes(void* A_data, ATimesFn ATimes);

    int setPreconditioner(void* P_data, PSetupFn Pset, PSolveFn Psol);

    int setScalingVectors(N_Vector s, N_Vector nul);

    int getNumIters();

    realtype getResNorm();

    N_Vector getResid();

};


class SUNLinSolSPBCGS : public SUNLinSolWrapper {
public:
    SUNLinSolSPBCGS(N_Vector y, int pretype, int maxl);


    SUNLinSolSPBCGS(AmiVector const& x, int pretype, int maxl);


    int setATimes(void* A_data, ATimesFn ATimes);

    int setPreconditioner(void* P_data, PSetupFn Pset, PSolveFn Psol);

    int setScalingVectors(N_Vector s, N_Vector nul);

    int getNumIters();

    realtype getResNorm();

    N_Vector getResid();

};


class SUNLinSolSPFGMR: public SUNLinSolWrapper {
public:
    /**
     * @brief SUNLinSolSPFGMR
     * @param y
     * @param pretype
     * @param maxl
     */
    SUNLinSolSPFGMR(AmiVector const& x, int pretype, int maxl);

    int setATimes(void* A_data, ATimesFn ATimes);

    int setPreconditioner(void* P_data, PSetupFn Pset, PSolveFn Psol);

    int setScalingVectors(N_Vector s, N_Vector nul);

    int getNumIters();

    realtype getResNorm();

    N_Vector getResid();

};


class SUNLinSolSPGMR: public SUNLinSolWrapper {
public:
    /**
     * @brief SUNLinSolSPGMR
     * @param y
     * @param pretype
     * @param maxl
     */
    SUNLinSolSPGMR(AmiVector const& x, int pretype, int maxl);


    int setATimes(void* A_data, ATimesFn ATimes);

    int setPreconditioner(void* P_data, PSetupFn Pset, PSolveFn Psol);

    int setScalingVectors(N_Vector s, N_Vector nul);

    int getNumIters();

    realtype getResNorm();

    N_Vector getResid();
};


class SUNLinSolSPTFQMR: public SUNLinSolWrapper {
public:
    /**
     * @brief SUNLinSolSPTFQMR
     * @param y
     * @param pretype
     * @param maxl
     */
    SUNLinSolSPTFQMR(N_Vector y, int pretype, int maxl);


    SUNLinSolSPTFQMR(AmiVector const& x, int pretype, int maxl);

    int setATimes(void* A_data, ATimesFn ATimes);

    int setPreconditioner(void* P_data, PSetupFn Pset, PSolveFn Psol);

    int setScalingVectors(N_Vector s, N_Vector nul);

    int getNumIters();

    realtype getResNorm();

    N_Vector getResid();

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

    int setSysFn(SUNNonlinSolSysFn SysFn);

    int setLSetupFn(SUNNonlinSolLSetupFn SetupFn);

    int setLSolveFn(SUNNonlinSolLSolveFn SolveFn);

    int setConvTestFn(SUNNonlinSolConvTestFn CTestFn);

    int setMaxIters(int maxiters);

    long int getNumIters();

    int getCurIter();

    long int getNumConvFails();

protected:

    /**
     * @brief initialize
     * @return
     */
    void initialize();

    /** the wrapper solver */
    SUNNonlinearSolver solver = nullptr;
};



class SUNNonLinSolNewton: public SUNNonLinSolWrapper {
public:

    SUNNonLinSolNewton(N_Vector x);


    SUNNonLinSolNewton(int count, N_Vector x);

    int getSysFn(SUNNonlinSolSysFn *SysFn);
};


class SUNNonLinSolFixedPoint: public SUNNonLinSolWrapper {
public:
    SUNNonLinSolFixedPoint(N_Vector x, int m = 0);

    SUNNonLinSolFixedPoint(int count, N_Vector x, int m = 0);

    int getSysFn(SUNNonlinSolSysFn *SysFn);
};


} // namespace amici
#endif // AMICI_SUNDIALS_LINSOL_WRAPPER_H

