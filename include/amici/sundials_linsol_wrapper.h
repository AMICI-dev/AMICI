#ifndef AMICI_SUNDIALS_LINSOL_WRAPPER_H
#define AMICI_SUNDIALS_LINSOL_WRAPPER_H

#include "amici/sundials_matrix_wrapper.h"
#include "amici/vector.h"

#include <sundials/sundials_config.h>
#include <sunlinsol/sunlinsol_band.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_klu.h>
#include <sunlinsol/sunlinsol_pcg.h>
#include <sunlinsol/sunlinsol_spbcgs.h>
#include <sunlinsol/sunlinsol_spfgmr.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunlinsol/sunlinsol_sptfqmr.h>
#ifdef SUNDIALS_SUPERLUMT
#include <sunlinsol/sunlinsol_superlumt.h>
#endif
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

namespace amici {

/**
 * @brief A RAII wrapper for SUNLinearSolver structs.
 *
 * For details on member functions see documentation in
 * sunlinsol/sundials_linearsolver.h.
 */
class SUNLinSolWrapper {
  public:
    SUNLinSolWrapper() = default;

    /**
     * @brief Wrap existing SUNLinearSolver
     * @param linsol
     */
    explicit SUNLinSolWrapper(SUNLinearSolver linsol);

    virtual ~SUNLinSolWrapper();

    /**
     * @brief Copy constructor
     * @param other
     */
    SUNLinSolWrapper(SUNLinSolWrapper const& other) = delete;

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
    SUNLinSolWrapper& operator=(SUNLinSolWrapper const& other) = delete;

    /**
     * @brief Move assignment
     * @param other
     * @return
     */
    SUNLinSolWrapper& operator=(SUNLinSolWrapper&& other) noexcept;

    /**
     * @brief Returns the wrapped SUNLinSol.
     * @return SUNLinearSolver
     */
    SUNLinearSolver get() const;

    /**
     * @brief Returns an identifier for the linear solver type.
     * @return
     */
    SUNLinearSolver_Type getType() const;

    /**
     * @brief Performs any linear solver setup needed, based on an updated
     * system matrix A.
     * @param A
     */
    void setup(SUNMatrix A) const;

    /**
     * @brief Performs any linear solver setup needed, based on an updated
     * system matrix A.
     * @param A
     */
    void setup(SUNMatrixWrapper const& A) const;

    /**
     * @brief Solves a linear system A*x = b
     * @param A
     * @param x A template for cloning vectors needed within the solver.
     * @param b
     * @param tol Tolerance (weighted 2-norm), iterative solvers only
     * @return error flag
     */
    int Solve(SUNMatrix A, N_Vector x, N_Vector b, realtype tol) const;

    /**
     * @brief Returns the last error flag encountered within the linear solver
     * @return error flag
     */
    long int getLastFlag() const;

    /**
     * @brief Returns the integer and real workspace sizes for the linear solver
     * @param lenrwLS output argument for size of real workspace
     * @param leniwLS output argument for size of integer workspace
     * @return workspace size
     */
    int space(long int* lenrwLS, long int* leniwLS) const;

    /**
     * @brief Get the matrix A (matrix solvers only).
     * @return A
     */
    virtual SUNMatrix getMatrix() const;

  protected:
    /**
     * @brief Performs linear solver initialization (assumes that all
     * solver-specific options have been set).
     * @return error code
     */
    int initialize();

    /** Wrapped solver */
    SUNLinearSolver solver_{nullptr};
};

/**
 * @brief SUNDIALS band direct solver.
 */
class SUNLinSolBand : public SUNLinSolWrapper {
  public:
    /**
     * @brief Create solver using existing matrix A without taking ownership of
     * A.
     * @param x A template for cloning vectors needed within the solver.
     * @param A square matrix
     */
    SUNLinSolBand(N_Vector x, SUNMatrix A);

    /**
     * @brief Create new band solver and matrix A.
     * @param x A template for cloning vectors needed within the solver.
     * @param ubw upper bandwidth of band matrix A
     * @param lbw lower bandwidth of band matrix A
     */
    SUNLinSolBand(AmiVector const& x, int ubw, int lbw);

    SUNMatrix getMatrix() const override;

  private:
    /** Matrix A for solver, only if created by here. */
    SUNMatrixWrapper A_;
};

/**
 * @brief SUNDIALS dense direct solver.
 */
class SUNLinSolDense : public SUNLinSolWrapper {
  public:
    /**
     * @brief Create dense solver
     * @param x A template for cloning vectors needed within the solver.
     */
    explicit SUNLinSolDense(AmiVector const& x);

    SUNMatrix getMatrix() const override;

  private:
    /** Matrix A for solver, only if created by here. */
    SUNMatrixWrapper A_;
};

/**
 * @brief SUNDIALS KLU sparse direct solver.
 */
class SUNLinSolKLU : public SUNLinSolWrapper {
  public:
    /** KLU state reordering (different from SuperLUMT ordering!) */
    enum class StateOrdering { AMD, COLAMD, natural };

    /**
     * @brief Create KLU solver with given matrix
     * @param x A template for cloning vectors needed within the solver.
     * @param A sparse matrix
     */
    SUNLinSolKLU(N_Vector x, SUNMatrix A);

    /**
     * @brief Create KLU solver and matrix to operate on
     * @param x A template for cloning vectors needed within the solver.
     * @param nnz Number of non-zeros in matrix A
     * @param sparsetype Sparse matrix type (CSC_MAT, CSR_MAT)
     * @param ordering
     */
    SUNLinSolKLU(
        AmiVector const& x, int nnz, int sparsetype, StateOrdering ordering
    );

    SUNMatrix getMatrix() const override;

    /**
     * @brief Reinitializes memory and flags for a new factorization
     * (symbolic and numeric) to be conducted at the next solver setup call.
     *
     * For more details see sunlinsol/sunlinsol_klu.h
     * @param nnz Number of non-zeros
     * @param reinit_type SUNKLU_REINIT_FULL or SUNKLU_REINIT_PARTIAL
     */
    void reInit(int nnz, int reinit_type);

    /**
     * @brief Sets the ordering used by KLU for reducing fill in the linear
     * solve.
     * @param ordering
     */
    void setOrdering(StateOrdering ordering);

  private:
    /** Sparse matrix A for solver, only if created by here. */
    SUNMatrixWrapper A_;
};

#ifdef SUNDIALS_SUPERLUMT
/**
 * @brief SUNDIALS SuperLUMT sparse direct solver.
 */
class SUNLinSolSuperLUMT : public SUNLinSolWrapper {
  public:
    /** SuperLUMT ordering (different from KLU ordering!) */
    enum class StateOrdering {
        natural,
        minDegATA,
        minDegATPlusA,
        COLAMD,
    };

    /**
     * @brief Create SuperLUMT solver with given matrix
     * @param x A template for cloning vectors needed within the solver.
     * @param A sparse matrix
     * @param numThreads Number of threads to be used by SuperLUMT
     */
    SUNLinSolSuperLUMT(N_Vector x, SUNMatrix A, int numThreads);

    /**
     * @brief Create SuperLUMT solver and matrix to operate on
     *
     * Will set number of threads according to environment variable
     * AMICI_SUPERLUMT_NUM_THREADS. Will default to 1 thread if unset.
     *
     * @param x A template for cloning vectors needed within the solver.
     * @param nnz Number of non-zeros in matrix A
     * @param sparsetype Sparse matrix type (CSC_MAT, CSR_MAT)
     * @param ordering
     */
    SUNLinSolSuperLUMT(
        AmiVector const& x, int nnz, int sparsetype, StateOrdering ordering
    );

    /**
     * @brief Create SuperLUMT solver and matrix to operate on
     * @param x A template for cloning vectors needed within the solver.
     * @param nnz Number of non-zeros in matrix A
     * @param sparsetype Sparse matrix type (CSC_MAT, CSR_MAT)
     * @param ordering
     * @param numThreads Number of threads to be used by SuperLUMT
     */
    SUNLinSolSuperLUMT(
        AmiVector const& x, int nnz, int sparsetype, StateOrdering ordering,
        int numThreads
    );

    SUNMatrix getMatrix() const override;

    /**
     * @brief Sets the ordering used by SuperLUMT for reducing fill in the
     * linear solve.
     * @param ordering
     */
    void setOrdering(StateOrdering ordering);

  private:
    /** Sparse matrix A for solver, only if created by here. */
    SUNMatrixWrapper A;
};

#endif

/**
 * @brief SUNDIALS scaled preconditioned CG (Conjugate Gradient method) (PCG)
 * solver.
 */
class SUNLinSolPCG : public SUNLinSolWrapper {
  public:
    /**
     * @brief Create PCG solver.
     * @param y
     * @param pretype Preconditioner type (PREC_NONE, PREC_LEFT, PREC_RIGHT,
     * PREC_BOTH)
     * @param maxl Maximum number of solver iterations
     */
    SUNLinSolPCG(N_Vector y, int pretype, int maxl);

    /**
     * @brief Sets the function pointer for ATimes
     * (see sundials/sundials_linearsolver.h).
     * @param A_data
     * @param ATimes
     * @return
     */
    int setATimes(void* A_data, ATimesFn ATimes);

    /**
     * @brief Sets function pointers for PSetup and PSolve routines inside
     * of iterative linear solver objects
     * (see sundials/sundials_linearsolver.h).
     * @param P_data
     * @param Pset
     * @param Psol
     * @return
     */
    int setPreconditioner(void* P_data, PSetupFn Pset, PSolveFn Psol);

    /**
     * @brief Sets pointers to left/right scaling vectors for the linear
     * system solve (see sundials/sundials_linearsolver.h).
     * @param s
     * @param nul
     * @return
     */
    int setScalingVectors(N_Vector s, N_Vector nul);

    /**
     * @brief Returns the number of linear iterations performed in the last
     * 'Solve' call
     * @return Number of iterations
     */
    int getNumIters() const;

    /**
     * @brief Returns the final residual norm from the last 'Solve' call.
     * @return residual norm
     */
    realtype getResNorm() const;

    /**
     * @brief Get preconditioned initial residual
     * (see sundials/sundials_linearsolver.h).
     * @return
     */
    N_Vector getResid() const;
};

/**
 * @brief SUNDIALS scaled preconditioned Bi-CGStab (Bi-Conjugate Gradient
 * Stable method) (SPBCGS) solver.
 */
class SUNLinSolSPBCGS : public SUNLinSolWrapper {
  public:
    /**
     * @brief SUNLinSolSPBCGS
     * @param x A template for cloning vectors needed within the solver.
     * @param pretype Preconditioner type (PREC_NONE, PREC_LEFT, PREC_RIGHT,
     * PREC_BOTH)
     * @param maxl Maximum number of solver iterations
     */
    explicit SUNLinSolSPBCGS(
        N_Vector x, int pretype = PREC_NONE, int maxl = SUNSPBCGS_MAXL_DEFAULT
    );

    /**
     * @brief SUNLinSolSPBCGS
     * @param x A template for cloning vectors needed within the solver.
     * @param pretype Preconditioner type (PREC_NONE, PREC_LEFT, PREC_RIGHT,
     * PREC_BOTH)
     * @param maxl Maximum number of solver iterations
     */
    explicit SUNLinSolSPBCGS(
        AmiVector const& x, int pretype = PREC_NONE,
        int maxl = SUNSPBCGS_MAXL_DEFAULT
    );

    /**
     * @brief Sets the function pointer for ATimes
     * (see sundials/sundials_linearsolver.h).
     * @param A_data
     * @param ATimes
     * @return
     */
    int setATimes(void* A_data, ATimesFn ATimes);

    /**
     * @brief Sets function pointers for PSetup and PSolve routines inside
     * of iterative linear solver objects
     * (see sundials/sundials_linearsolver.h).
     * @param P_data
     * @param Pset
     * @param Psol
     * @return
     */
    int setPreconditioner(void* P_data, PSetupFn Pset, PSolveFn Psol);

    /**
     * @brief Sets pointers to left/right scaling vectors for the linear
     * system solve (see sundials/sundials_linearsolver.h).
     * @param s
     * @param nul
     * @return
     */
    int setScalingVectors(N_Vector s, N_Vector nul);

    /**
     * @brief Returns the number of linear iterations performed in the last
     * 'Solve' call
     * @return Number of iterations
     */
    int getNumIters() const;

    /**
     * @brief Returns the final residual norm from the last 'Solve' call.
     * @return residual norm
     */
    realtype getResNorm() const;

    /**
     * @brief Get preconditioned initial residual
     * (see sundials/sundials_linearsolver.h).
     * @return
     */
    N_Vector getResid() const;
};

/**
 * @brief SUNDIALS scaled preconditioned FGMRES (Flexible Generalized Minimal
 * Residual method) (SPFGMR) solver.
 */
class SUNLinSolSPFGMR : public SUNLinSolWrapper {
  public:
    /**
     * @brief SUNLinSolSPFGMR
     * @param x A template for cloning vectors needed within the solver.
     * @param pretype Preconditioner type (PREC_NONE, PREC_LEFT, PREC_RIGHT,
     * PREC_BOTH)
     * @param maxl Maximum number of solver iterations
     */
    SUNLinSolSPFGMR(AmiVector const& x, int pretype, int maxl);

    /**
     * @brief Sets the function pointer for ATimes
     * (see sundials/sundials_linearsolver.h).
     * @param A_data
     * @param ATimes
     * @return
     */
    int setATimes(void* A_data, ATimesFn ATimes);

    /**
     * @brief Sets function pointers for PSetup and PSolve routines inside
     * of iterative linear solver objects
     * (see sundials/sundials_linearsolver.h).
     * @param P_data
     * @param Pset
     * @param Psol
     * @return
     */
    int setPreconditioner(void* P_data, PSetupFn Pset, PSolveFn Psol);

    /**
     * @brief Sets pointers to left/right scaling vectors for the linear
     * system solve (see sundials/sundials_linearsolver.h).
     * @param s
     * @param nul
     * @return
     */
    int setScalingVectors(N_Vector s, N_Vector nul);

    /**
     * @brief Returns the number of linear iterations performed in the last
     * 'Solve' call
     * @return Number of iterations
     */
    int getNumIters() const;

    /**
     * @brief Returns the final residual norm from the last 'Solve' call.
     * @return residual norm
     */
    realtype getResNorm() const;

    /**
     * @brief Get preconditioned initial residual
     * (see sundials/sundials_linearsolver.h).
     * @return
     */
    N_Vector getResid() const;
};

/**
 * @brief SUNDIALS scaled preconditioned GMRES (Generalized Minimal Residual
 * method) solver (SPGMR).
 */
class SUNLinSolSPGMR : public SUNLinSolWrapper {
  public:
    /**
     * @brief Create SPGMR solver
     * @param x A template for cloning vectors needed within the solver.
     * @param pretype Preconditioner type (PREC_NONE, PREC_LEFT, PREC_RIGHT,
     * PREC_BOTH)
     * @param maxl Maximum number of solver iterations
     */
    explicit SUNLinSolSPGMR(
        AmiVector const& x, int pretype = PREC_NONE,
        int maxl = SUNSPGMR_MAXL_DEFAULT
    );

    /**
     * @brief Sets the function pointer for ATimes
     * (see sundials/sundials_linearsolver.h).
     * @param A_data
     * @param ATimes
     * @return
     */
    int setATimes(void* A_data, ATimesFn ATimes);

    /**
     * @brief Sets function pointers for PSetup and PSolve routines inside
     * of iterative linear solver objects
     * (see sundials/sundials_linearsolver.h).
     * @param P_data
     * @param Pset
     * @param Psol
     * @return
     */
    int setPreconditioner(void* P_data, PSetupFn Pset, PSolveFn Psol);

    /**
     * @brief Sets pointers to left/right scaling vectors for the linear
     * system solve (see sundials/sundials_linearsolver.h).
     * @param s
     * @param nul
     * @return
     */
    int setScalingVectors(N_Vector s, N_Vector nul);

    /**
     * @brief Returns the number of linear iterations performed in the last
     * 'Solve' call
     * @return Number of iterations
     */
    int getNumIters() const;

    /**
     * @brief Returns the final residual norm from the last 'Solve' call.
     * @return residual norm
     */
    realtype getResNorm() const;

    /**
     * @brief Get preconditioned initial residual
     * (see sundials/sundials_linearsolver.h).
     * @return
     */
    N_Vector getResid() const;
};

/**
 * @brief SUNDIALS scaled preconditioned TFQMR (Transpose-Free Quasi-Minimal
 * Residual method) (SPTFQMR) solver.
 */
class SUNLinSolSPTFQMR : public SUNLinSolWrapper {
  public:
    /**
     * @brief Create SPTFQMR solver
     * @param x A template for cloning vectors needed within the solver.
     * @param pretype Preconditioner type (PREC_NONE, PREC_LEFT, PREC_RIGHT,
     * PREC_BOTH)
     * @param maxl Maximum number of solver iterations
     */
    explicit SUNLinSolSPTFQMR(
        N_Vector x, int pretype = PREC_NONE, int maxl = SUNSPTFQMR_MAXL_DEFAULT
    );

    /**
     * @brief Create SPTFQMR solver
     * @param x A template for cloning vectors needed within the solver.
     * @param pretype Preconditioner type (PREC_NONE, PREC_LEFT, PREC_RIGHT,
     * PREC_BOTH)
     * @param maxl Maximum number of solver iterations
     */
    explicit SUNLinSolSPTFQMR(
        AmiVector const& x, int pretype = PREC_NONE,
        int maxl = SUNSPTFQMR_MAXL_DEFAULT
    );

    /**
     * @brief Sets the function pointer for ATimes
     * (see sundials/sundials_linearsolver.h).
     * @param A_data
     * @param ATimes
     * @return
     */
    int setATimes(void* A_data, ATimesFn ATimes);

    /**
     * @brief Sets function pointers for PSetup and PSolve routines inside
     * of iterative linear solver objects
     * (see sundials/sundials_linearsolver.h).
     * @param P_data
     * @param Pset
     * @param Psol
     * @return
     */
    int setPreconditioner(void* P_data, PSetupFn Pset, PSolveFn Psol);

    /**
     * @brief Sets pointers to left/right scaling vectors for the linear
     * system solve (see sundials/sundials_linearsolver.h).
     * @param s
     * @param nul
     * @return
     */
    int setScalingVectors(N_Vector s, N_Vector nul);

    /**
     * @brief Returns the number of linear iterations performed in the last
     * 'Solve' call
     * @return Number of iterations
     */
    int getNumIters() const;

    /**
     * @brief Returns the final residual norm from the last 'Solve' call.
     * @return residual norm
     */
    realtype getResNorm() const;

    /**
     * @brief Get preconditioned initial residual
     * (see sundials/sundials_linearsolver.h).
     * @return
     */
    N_Vector getResid() const;
};

/**
 * @brief A RAII wrapper for SUNNonLinearSolver structs which solve the
 * nonlinear system F (y) = 0 or G(y) = y.
 */
class SUNNonLinSolWrapper {
  public:
    /**
     * @brief SUNNonLinSolWrapper from existing SUNNonlinearSolver
     * @param sol
     */
    explicit SUNNonLinSolWrapper(SUNNonlinearSolver sol);

    virtual ~SUNNonLinSolWrapper();

    /**
     * @brief Copy constructor
     * @param other
     */
    SUNNonLinSolWrapper(SUNNonLinSolWrapper const& other) = delete;

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
    SUNNonLinSolWrapper& operator=(SUNNonLinSolWrapper const& other) = delete;

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
     * @brief Get type ID of the solver
     * @return
     */
    SUNNonlinearSolver_Type getType() const;

    /**
     * @brief Setup solver
     * @param y  the initial iteration passed to the nonlinear solver.
     * @param mem the sundials integrator memory structure.
     * @return
     */
    int setup(N_Vector y, void* mem);

    /**
     * @brief Solve the nonlinear system F (y) = 0 or G(y) = y.
     * @param y0 the initial iterate for the nonlinear solve. This must remain
     * unchanged throughout the solution process.
     * @param y the solution to the nonlinear system
     * @param w the solution error weight vector used for computing weighted
     * error norms.
     * @param tol the requested solution tolerance in the weighted root-mean-
     * squared norm.
     * @param callLSetup a flag indicating that the integrator recommends for
     * the linear solver setup function to be called.
     * @param mem the sundials integrator memory structure.
     * @return
     */
    int Solve(
        N_Vector y0, N_Vector y, N_Vector w, realtype tol, bool callLSetup,
        void* mem
    );

    /**
     * @brief Set function to evaluate the nonlinear residual function F(y) = 0
     * or the fixed point function G(y) = y
     * @param SysFn
     * @return
     */
    int setSysFn(SUNNonlinSolSysFn SysFn);

    /**
     * @brief Set linear solver setup function.
     * @param SetupFn
     * @return
     */
    int setLSetupFn(SUNNonlinSolLSetupFn SetupFn);

    /**
     * @brief Set linear solver solve function.
     * @param SolveFn
     * @return
     */
    int setLSolveFn(SUNNonlinSolLSolveFn SolveFn);

    /**
     * @brief Set function to test for convergence
     * @param CTestFn
     * @param ctest_data
     * @return
     */
    int setConvTestFn(SUNNonlinSolConvTestFn CTestFn, void* ctest_data);

    /**
     * @brief Set maximum number of non-linear iterations
     * @param maxiters
     * @return
     */
    int setMaxIters(int maxiters);

    /**
     * @brief getNumIters
     * @return
     */
    long int getNumIters() const;

    /**
     * @brief getCurIter
     * @return
     */
    int getCurIter() const;

    /**
     * @brief getNumConvFails
     * @return
     */
    long int getNumConvFails() const;

  protected:
    /**
     * @brief initialize
     */
    void initialize();

    /** the wrapper solver */
    SUNNonlinearSolver solver = nullptr;
};

/**
 * @brief SUNDIALS Newton non-linear solver to solve F (y) = 0.
 */
class SUNNonLinSolNewton : public SUNNonLinSolWrapper {
  public:
    /**
     * @brief Create Newton solver
     * @param x A template for cloning vectors needed within the solver.
     */
    explicit SUNNonLinSolNewton(N_Vector x);

    /**
     * @brief Create Newton solver for enabled sensitivity analysis
     * @param count Number of vectors in the nonlinear solve. When integrating
     * a system containing Ns sensitivities the value of count is:
     *    - Ns+1 if using a simultaneous corrector approach.
     *    - Ns if using a staggered corrector approach.
     * @param x A template for cloning vectors needed within the solver.
     */
    SUNNonLinSolNewton(int count, N_Vector x);

    /**
     * @brief Get function to evaluate the nonlinear residual function F(y) = 0
     * @param SysFn
     * @return
     */
    int getSysFn(SUNNonlinSolSysFn* SysFn) const;
};

/**
 * @brief SUNDIALS Fixed point non-linear solver to solve G(y) = y.
 */
class SUNNonLinSolFixedPoint : public SUNNonLinSolWrapper {
  public:
    /**
     * @brief Create fixed-point solver
     * @param x template for cloning vectors needed within the solver.
     * @param m number of acceleration vectors to use
     */
    explicit SUNNonLinSolFixedPoint(const_N_Vector x, int m = 0);

    /**
     * @brief Create fixed-point solver for use with sensitivity analysis
     * @param count Number of vectors in the nonlinear solve. When integrating
     * a system containing Ns sensitivities the value of count is:
     *    - Ns+1 if using a simultaneous corrector approach.
     *    - Ns if using a staggered corrector approach.
     * @param x template for cloning vectors needed within the solver.
     * @param m number of acceleration vectors to use
     */
    SUNNonLinSolFixedPoint(int count, const_N_Vector x, int m = 0);

    /**
     * @brief Get function to evaluate the fixed point function G(y) = y
     * @param SysFn
     * @return
     */
    int getSysFn(SUNNonlinSolSysFn* SysFn) const;
};

} // namespace amici
#endif // AMICI_SUNDIALS_LINSOL_WRAPPER_H
