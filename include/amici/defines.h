#ifndef AMICI_DEFINES_H
#define AMICI_DEFINES_H
#include <cmath>

/** MS definition of PI and other constants */
#define _USE_MATH_DEFINES
#ifndef M_PI
/** define PI if we still have no definition */
#define M_PI 3.14159265358979323846
#endif

namespace amici {

constexpr double pi = M_PI;

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

typedef double realtype;

typedef enum {
    AMICI_BLAS_RowMajor = 101,
    AMICI_BLAS_ColMajor = 102
} AMICI_BLAS_LAYOUT;

typedef enum {
    AMICI_BLAS_NoTrans = 111,
    AMICI_BLAS_Trans = 112,
    AMICI_BLAS_ConjTrans = 113
} AMICI_BLAS_TRANSPOSE;

typedef enum AMICI_parameter_scaling_TAG {
    AMICI_SCALING_NONE,
    AMICI_SCALING_LN,
    AMICI_SCALING_LOG10
} AMICI_parameter_scaling;

typedef enum AMICI_o2mode_TAG {
    AMICI_O2MODE_NONE,
    AMICI_O2MODE_FULL,
    AMICI_O2MODE_DIR
} AMICI_o2mode;

typedef enum AMICI_sensi_order_TAG {
    AMICI_SENSI_ORDER_NONE,
    AMICI_SENSI_ORDER_FIRST,
    AMICI_SENSI_ORDER_SECOND
} AMICI_sensi_order;

typedef enum AMICI_sensi_meth_TAG {
    AMICI_SENSI_NONE,
    AMICI_SENSI_FSA,
    AMICI_SENSI_ASA,
    AMICI_SENSI_SS
} AMICI_sensi_meth;

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


enum InternalSensitivityMethod {
    SIMULTANEOUS = 1,
    STAGGERED = 2,
    STAGGERED1 = 3
};

enum InterpolationType {
    HERMITE = 1,
    POLYNOMIAL = 2
};

enum LinearMultistepMethod {
    ADAMS = 1,
    BDF = 2
};

enum NonlinearSolverIteration {
    FUNCTIONAL = 1,
    NEWTON = 2
};

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
