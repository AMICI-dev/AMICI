#ifndef amici_h
#define amici_h

#include <include/symbolic_functions.h>
#include <include/udata.h>
#include <include/rdata.h>
#include <include/edata.h>
#include <include/tdata.h>
#include <include/amici_solver.h>
#include <cstdbool>
#include <cvodes/cvodes.h>

#include <include/amici_defines.h>

// ensure definitions are in sync
static_assert(AMICI_SUCCESS == CV_SUCCESS, "AMICI_SUCCESS != CV_SUCCESS");
static_assert(AMICI_DATA_RETURN == CV_TSTOP_RETURN, "AMICI_DATA_RETURN != CV_TSTOP_RETURN");
static_assert(AMICI_ROOT_RETURN == CV_ROOT_RETURN, "AMICI_ROOT_RETURN != CV_ROOT_RETURN");
static_assert(AMICI_NORMAL == CV_NORMAL, "AMICI_NORMAL != CV_NORMAL");
static_assert(AMICI_ONE_STEP == CV_ONE_STEP, "AMICI_ONE_STEP != CV_ONE_STEP");

void printErrMsgIdAndTxt(const char * identifier, const char *msg, ...);

void printWarnMsgIdAndTxt(const char * identifier, const char *msg, ...);

/**
 * @brief msgIdAndTxtFp
 * @param identifier string with error message identifier
 * @param err_msg string with error message printf-style format
 * @param ... unused
 */
typedef void (*msgIdAndTxtFp)(const char * identifier, const char * err_msg, ...);

// function pointers to process errors / warnings
extern msgIdAndTxtFp errMsgIdAndTxt;
extern msgIdAndTxtFp warnMsgIdAndTxt;

int runAmiciSimulation(UserData *udata, const ExpData *edata, ReturnData *rdata, Model *model, Solver *solver);

int prepDataSensis(int it, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata);

int getDataSensisFSA(int it, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata, Solver *solver, Model *model);

int getDataOutput(int it, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata, Solver *solver, Model *model);


realtype getTnext(realtype *troot, int iroot, realtype *tdata, int it, Model *model);


void amici_dgemv(AMICI_BLAS_LAYOUT layout,
                 AMICI_BLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY);

void amici_dgemm(AMICI_BLAS_LAYOUT layout, AMICI_BLAS_TRANSPOSE TransA,
                 AMICI_BLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);

#endif /* amici_h */
