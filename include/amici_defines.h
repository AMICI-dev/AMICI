#ifndef AMICI_DEFINES_H
#define AMICI_DEFINES_H

// clang-format off

#define AMICI_ONEOUTPUT   5

/* Return codes */
#define AMICI_ERROR_UDATA              -99
#define AMICI_ERROR_EDATA              -98
#define AMICI_ERROR_RDATA              -97
#define AMICI_ERROR_TDATA              -96
#define AMICI_ERROR_SETUP              -95
#define AMICI_ERROR_SETUPB             -94
#define AMICI_ERROR_NOTHINGTODO        -93
#define AMICI_ERROR_FSA                -92
#define AMICI_ERROR_ASA                -91
#define AMICI_ERROR_SA                 -90
#define AMICI_ERROR_SS_SENSIS          -89
#define AMICI_ERROR_DATA               -88
#define AMICI_ERROR_EVENT              -87
#define AMICI_ERROR_SIMULATION         -86
#define AMICI_ERROR_NEWTONSOLVER       -85
#define AMICI_ERROR_NEWTON_LINSOLVER   -84
#define AMICI_ERROR_NOT_IMPLEMENTED    -83
#define AMICI_ERROR_MODEL              -82
#define AMICI_ERROR_OTHER              -81
#define AMICI_ERROR_SIM2STEADYSTATE    -80
#define AMICI_ERROR_PREEQUILIBRATION   -79
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
// clang-format on

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


#endif
