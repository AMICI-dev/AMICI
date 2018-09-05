#ifndef amici_h
#define amici_h

#include "amici/model.h"
#include "amici/solver.h"
#include "amici/exception.h"
#include "amici/defines.h"
#include "amici/rdata.h"
#include "amici/edata.h"
#include "amici/symbolic_functions.h"
#include "amici/cblas.h"

namespace amici {

void printErrMsgIdAndTxt(const char *identifier, const char *format, ...);

void printWarnMsgIdAndTxt(const char *identifier, const char *format, ...);

// function pointers to process errors / warnings
extern msgIdAndTxtFp errMsgIdAndTxt;
extern msgIdAndTxtFp warnMsgIdAndTxt;

std::unique_ptr<ReturnData> runAmiciSimulation(Solver &solver, const ExpData *edata, Model &model);

} // namespace amici

#endif /* amici_h */
