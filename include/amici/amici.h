#ifndef amici_h
#define amici_h

#include "amici/cblas.h"
#include "amici/defines.h"
#include "amici/edata.h"
#include "amici/exception.h"
#include "amici/model.h"
#include "amici/rdata.h"
#include "amici/solver.h"
#include "amici/symbolic_functions.h"

namespace amici {

/*!
 * @brief printErrMsgIdAndTxt prints a specified error message associated to the
 * specified identifier
 *
 * @param[in] identifier error identifier @type char
 * @param[in] format string with error message printf-style format
 * @param ... arguments to be formatted
 * @return void
 */
void
printErrMsgIdAndTxt(const char* identifier, const char* format, ...);

/*!
 * @brief printErrMsgIdAndTxt prints a specified warning message associated to
 * the specified identifier
 *
 * @param[in] identifier warning identifier @type char
 * @param[in] format string with error message printf-style format
 * @param ... arguments to be formatted
 * @return void
 */
void
printWarnMsgIdAndTxt(const char* identifier, const char* format, ...);

/** errMsgIdAndTxt is a function pointer for printErrMsgIdAndTxt  */
extern msgIdAndTxtFp errMsgIdAndTxt;

/** warnMsgIdAndTxt is a function pointer for printWarnMsgIdAndTxt  */
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
std::unique_ptr<ReturnData>
runAmiciSimulation(Solver& solver,
                   const ExpData* edata,
                   Model& model,
                   bool rethrow = false);

/*!
 * runAmiciSimulations does the same as runAmiciSimulation, but for multiple
 * ExpData instances.
 *
 * @param solver Solver instance
 * @param edatas experimental data objects
 * @param model model specification object
 * @param failfast flag to allow early termination
 * @param num_threads number of threads for parallel execution
 * @return vector of pointers to return data objects
 */
std::vector<std::unique_ptr<ReturnData>>
runAmiciSimulations(Solver const& solver,
                    const std::vector<ExpData*>& edatas,
                    Model const& model,
                    bool failfast,
                    int num_threads);

} // namespace amici

#endif /* amici_h */
