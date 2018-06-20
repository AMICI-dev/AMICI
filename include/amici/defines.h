#ifndef AMICI_DEFINES_H
#define AMICI_DEFINES_H
#include <cmath>

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
#define AMICI_TOO_MUCH_WORK             -1
#define AMICI_TOO_MUCH_ACC              -2
#define AMICI_ERR_FAILURE               -3
#define AMICI_CONV_FAILURE              -4
#define AMICI_ILL_INPUT                -22
#define AMICI_ERROR                    -99
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

/** defines variable type for simulation variables (determines numerical accuracy) */
typedef double realtype;

/** BLAS Matrix Layout, affects dgemm and gemv calls */
typedef enum {
    AMICI_BLAS_RowMajor = 101,
    AMICI_BLAS_ColMajor = 102
} AMICI_BLAS_LAYOUT;

/** BLAS Matrix Transposition, affects dgemm and gemv calls */
typedef enum {
    AMICI_BLAS_NoTrans = 111,
    AMICI_BLAS_Trans = 112,
    AMICI_BLAS_ConjTrans = 113
} AMICI_BLAS_TRANSPOSE;

/** modes for parameter transformations */
typedef enum AMICI_parameter_scaling_TAG {
    AMICI_SCALING_NONE,
    AMICI_SCALING_LN,
    AMICI_SCALING_LOG10
} AMICI_parameter_scaling;

/** modes for second order sensitivity analysis */
typedef enum AMICI_o2mode_TAG {
    AMICI_O2MODE_NONE,
    AMICI_O2MODE_FULL,
    AMICI_O2MODE_DIR
} AMICI_o2mode;

/** orders of sensitivity analysis */
typedef enum AMICI_sensi_order_TAG {
    AMICI_SENSI_ORDER_NONE,
    AMICI_SENSI_ORDER_FIRST,
    AMICI_SENSI_ORDER_SECOND
} AMICI_sensi_order;

/** methods for sensitivity computation */
typedef enum AMICI_sensi_meth_TAG {
    AMICI_SENSI_NONE,
    AMICI_SENSI_FSA,
    AMICI_SENSI_ASA
} AMICI_sensi_meth;

/** linear solvers for CVODES/IDAS */
enum LinearSolver {
    AMICI_DENSE       = 1,
    AMICI_BAND        = 2,
    AMICI_LAPACKDENSE = 3,
    AMICI_LAPACKBAND  = 4,
    AMICI_DIAG        = 5,
    AMICI_SPGMR       = 6,
    AMICI_SPBCG       = 7,
    AMICI_SPTFQMR     = 8,
    AMICI_KLU         = 9
};

/** CVODES/IDAS forward sensitivity computation method */
enum InternalSensitivityMethod {
    SIMULTANEOUS = 1,
    STAGGERED = 2,
    STAGGERED1 = 3
};

/** CVODES/IDAS state interpolation for adjoint sensitivity analysis */
enum InterpolationType {
    HERMITE = 1,
    POLYNOMIAL = 2
};

/** CVODES/IDAS linear multistep method */
enum LinearMultistepMethod {
    ADAMS = 1,
    BDF = 2
};

/** CVODES/IDAS Nonlinear Iteration method */
enum NonlinearSolverIteration {
    FUNCTIONAL = 1,
    NEWTON = 2
};

/** KLU state reordering */
enum StateOrdering {
    AMD,
    COLAMD,
    natural
};
    
    /**
     * @brief msgIdAndTxtFp
     * @param identifier string with error message identifier
     * @param format string with error message printf-style format
     * @param ... arguments to be formatted
     */
    typedef void (*msgIdAndTxtFp)(const char *identifier, const char *format, ...);

// clang-format on

} // namespace amici

#endif
