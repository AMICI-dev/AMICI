#ifndef amici_h
#define amici_h

#include <cvodes/cvodes.h>
#include <include/symbolic_functions.h>
#include <include/amici_defines.h>

namespace amici {

class UserData;
class ReturnData;
class ExpData;
class Model;

// ensure definitions are in sync
static_assert(AMICI_SUCCESS == CV_SUCCESS, "AMICI_SUCCESS != CV_SUCCESS");
static_assert(AMICI_DATA_RETURN == CV_TSTOP_RETURN,
              "AMICI_DATA_RETURN != CV_TSTOP_RETURN");
static_assert(AMICI_ROOT_RETURN == CV_ROOT_RETURN,
              "AMICI_ROOT_RETURN != CV_ROOT_RETURN");
static_assert(AMICI_ILL_INPUT == CV_ILL_INPUT,
              "AMICI_ILL_INPUT != CV_ILL_INPUT");
    
    
    
static_assert(AMICI_NORMAL == CV_NORMAL, "AMICI_NORMAL != CV_NORMAL");
static_assert(AMICI_ONE_STEP == CV_ONE_STEP, "AMICI_ONE_STEP != CV_ONE_STEP");


void printErrMsgIdAndTxt(const char *identifier, const char *format, ...);

void printWarnMsgIdAndTxt(const char *identifier, const char *format, ...);

// function pointers to process errors / warnings
extern msgIdAndTxtFp errMsgIdAndTxt;
extern msgIdAndTxtFp warnMsgIdAndTxt;


void runAmiciSimulation(const UserData *udata, const ExpData *edata,
                       ReturnData *rdata, Model *model);

void amici_dgemv(AMICI_BLAS_LAYOUT layout, AMICI_BLAS_TRANSPOSE TransA,
                 const int M, const int N, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY);

void amici_dgemm(AMICI_BLAS_LAYOUT layout, AMICI_BLAS_TRANSPOSE TransA,
                 AMICI_BLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);

} // namespace amici

#endif /* amici_h */
