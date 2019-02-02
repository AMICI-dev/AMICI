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

/*!
 * runAmiciSimulation is the core integration routine. It initializes the solver
 * and runs the forward and backward problem.
 *
 * @param solver Solver instance
 * @param edata pointer to experimental data object
 * @param model model specification object
 * @param rethrow rethrow integration exceptions?
 * @return rdata pointer to return data object
 */
std::unique_ptr<ReturnData> runAmiciSimulation(Solver &solver, const ExpData *edata, Model &model, bool rethrow=false);

/*!
 * runAmiciSimulations does the same as runAmiciSimulation, but for multiple ExpData instances.
 *
 * @param solver Solver instance
 * @param edatas experimental data objects
 * @param model model specification object
 * @param num_threads number of threads for parallel execution
 * @return vector of pointers to return data objects
 */
std::vector<std::unique_ptr<ReturnData>> runAmiciSimulations(Solver const& solver,
                                                             const std::vector<ExpData *> &edatas,
                                                             Model const& model,
                                                             int num_threads);

} // namespace amici

#endif /* amici_h */
