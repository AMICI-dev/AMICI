#ifndef AMICI_DEFINES_H
#define AMICI_DEFINES_H

#include <cmath>
#include <functional>
#include <string>

namespace amici {

#define _USE_MATH_DEFINES
#ifdef M_PI
/** pi definition from MATH_DEFINES */
constexpr double pi = M_PI;
#else
/** MS definition of PI and other constants */
constexpr double pi = 3.14159265358979323846;
#endif

// clang-format off

constexpr int AMICI_ONEOUTPUT=                 5;

/* Return codes */
constexpr int AMICI_RECOVERABLE_ERROR=         1;
constexpr int AMICI_UNRECOVERABLE_ERROR=     -10;
constexpr int AMICI_TOO_MUCH_WORK=            -1;
constexpr int AMICI_TOO_MUCH_ACC=             -2;
constexpr int AMICI_ERR_FAILURE=              -3;
constexpr int AMICI_CONV_FAILURE=             -4;
constexpr int AMICI_ILL_INPUT=               -22;
constexpr int AMICI_ERROR=                   -99;
constexpr int AMICI_DAMPING_FACTOR_ERROR=   -800;
constexpr int AMICI_SINGULAR_JACOBIAN=      -807;
constexpr int AMICI_NOT_IMPLEMENTED=        -999;
constexpr int AMICI_SUCCESS=                   0;
constexpr int AMICI_DATA_RETURN=               1;
constexpr int AMICI_ROOT_RETURN=               2;

constexpr int AMICI_NORMAL=                    1;
constexpr int AMICI_ONE_STEP=                  2;

constexpr int AMICI_PREEQUILIBRATE=           -1;

#ifndef booleantype
#define booleantype int
#endif

/** defines variable type for simulation variables
 * (determines numerical accuracy) */
using realtype = double;

/** BLAS Matrix Layout, affects dgemm and gemv calls */
enum class BLASLayout{
    rowMajor = 101,
    colMajor = 102
};

/** BLAS Matrix Transposition, affects dgemm and gemv calls */
enum class BLASTranspose {
    noTrans = 111,
    trans = 112,
    conjTrans = 113
};

/** modes for parameter transformations */
enum class ParameterScaling {
    none,
    ln,
    log10
};

/** modes for second order sensitivity analysis */
enum class SecondOrderMode {
    none,
    full,
    directional
};

/** orders of sensitivity analysis */
enum class SensitivityOrder {
    none,
    first,
    second
};

/** methods for sensitivity computation */
enum class SensitivityMethod {
    none,
    forward,
    adjoint
};

/** linear solvers for CVODES/IDAS */
enum class LinearSolver {
    dense       = 1,
    band        = 2,
    LAPACKDense = 3,
    LAPACKBand  = 4,
    diag        = 5,
    SPGMR       = 6,
    SPBCG       = 7,
    SPTFQMR     = 8,
    KLU         = 9,
    SuperLUMT   = 10,
};

/** CVODES/IDAS forward sensitivity computation method */
enum class InternalSensitivityMethod {
    simultaneous = 1,
    staggered = 2,
    staggered1 = 3
};

/** CVODES/IDAS state interpolation for adjoint sensitivity analysis */
enum class InterpolationType {
    hermite = 1,
    polynomial = 2
};

/** CVODES/IDAS linear multistep method */
enum class LinearMultistepMethod {
    adams = 1,
    BDF = 2
};

/** CVODES/IDAS Nonlinear Iteration method */
enum class NonlinearSolverIteration {
    functional = 1, /** deprecated */
    fixedpoint = 1,
    newton = 2
};

/** Sensitivity computation mode in steadyStateProblem */
enum class SteadyStateSensitivityMode {
    newtonOnly,
    simulationFSA
};

/** State in which the steady state computation finished */
enum class SteadyStateStatus {
    failed_damping = -4,
    failed_factorization = -3,
    failed_convergence = -2,
    failed = -1,
    not_run = 0,
    success = 1
};

/** Context for which the sensitivity flag should be computed */
enum class SteadyStateContext {
    newton = 0,
    storage = 1,
    integration = 2
};

/** Damping factor flag for the Newton method */
enum class NewtonDampingFactorMode {
    off = 0,
    on = 1
};

/** fixedParameter to be used in condition context */
enum class FixedParameterContext {
    simulation = 0,
    preequilibration = 1,
    presimulation = 2,
};

enum class RDataReporting {
    full,
    residuals,
    likelihood,
};

/**
 * Type for function to process warnings or error messages.
 */
using outputFunctionType = std::function<void(std::string const& identifier,
                                              std::string const& message)>;

// clang-format on

} // namespace amici

#endif
