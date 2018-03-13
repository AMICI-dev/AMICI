#ifndef amici_h
#define amici_h

#include <include/symbolic_functions.h>
#include <include/amici_defines.h>
#include <include/amici_model.h>
#include <include/amici_solver.h>
#include "include/amici_interface_cpp.h"

namespace amici {

class ReturnData;
class ExpData;
class Model;
class Solver;

void printErrMsgIdAndTxt(const char *identifier, const char *format, ...);

void printWarnMsgIdAndTxt(const char *identifier, const char *format, ...);

// function pointers to process errors / warnings
extern msgIdAndTxtFp errMsgIdAndTxt;
extern msgIdAndTxtFp warnMsgIdAndTxt;


std::unique_ptr<ReturnData> runAmiciSimulation(Solver &solver, const ExpData *edata, Model &model);

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
