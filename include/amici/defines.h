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

#define AMICI_ONEOUTPUT   5

/* Return codes */
#define AMICI_RECOVERABLE_ERROR          1
#define AMICI_UNRECOVERABLE_ERROR      -10
#define AMICI_TOO_MUCH_WORK             -1
#define AMICI_TOO_MUCH_ACC              -2
#define AMICI_ERR_FAILURE               -3
#define AMICI_CONV_FAILURE              -4
#define AMICI_ILL_INPUT                -22
#define AMICI_ERROR                    -99
#define AMICI_NOT_IMPLEMENTED         -999
#define AMICI_SUCCESS                    0
#define AMICI_DATA_RETURN                1
#define AMICI_ROOT_RETURN                2

#define AMICI_NORMAL                     1
#define AMICI_ONE_STEP                   2

#define AMICI_PREEQUILIBRATE            -1

#ifndef booleantype
#define booleantype int
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
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
enum class NewtonStatus {
    failed=-1,
    newt=1,
    newt_sim=2,
    newt_sim_newt=3,
};

/** Damping factor flag for the Newton method */
enum class NewtonDampingFactorMode {
    off = 0,
    on = 1
};

/**
 * Type for function to process warnings or error messages.
 */
using outputFunctionType = std::function<void(std::string const& identifier,
                                              std::string const& message)>;

// clang-format on

} // namespace amici

#endif
