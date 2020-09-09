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
 * @brief Prints a specified error message associated with the specified
 * identifier
 *
 * @param id error identifier
 * @param message error message
 */
void printErrMsgIdAndTxt(std::string const &id, std::string const &message);

/*!
 * @brief Prints a specified warning message associated with the specified
 * identifier
 *
 * @param id warning identifier
 * @param message warning message
 */
void printWarnMsgIdAndTxt(std::string const &id, std::string const &message);

/**
 * @brief Main class for making calls to AMICI.
 *
 * This class is used to provide separate AMICI contexts, for example, for use
 * in multi-threaded applications where different threads want to use AMICI with
 * different settings, such custom logging functions.
 *
 * NOTE: For this moment, the context object needs to be set manually to any
 * Model and Solver object. If not set, they will use the default output
 * channel.
 */
class AmiciApplication {
  public:
    AmiciApplication() = default;

    /**
     * @brief Core integration routine. Initializes the solver and runs the
     * forward and backward problem.
     *
     * @param solver Solver instance
     * @param edata pointer to experimental data object
     * @param model model specification object
     * @param rethrow rethrow integration exceptions?
     * @return rdata pointer to return data object
     */
    std::unique_ptr<ReturnData> runAmiciSimulation(Solver &solver,
                                                   const ExpData *edata,
                                                   Model &model,
                                                   bool rethrow = false);

    /**
     * @brief Same as runAmiciSimulation, but for multiple ExpData instances.
     *
     * @param solver Solver instance
     * @param edatas experimental data objects
     * @param model model specification object
     * @param failfast flag to allow early termination
     * @param num_threads number of threads for parallel execution
     * @return vector of pointers to return data objects
     */
    std::vector<std::unique_ptr<ReturnData>>
    runAmiciSimulations(Solver const &solver,
                        const std::vector<ExpData *> &edatas,
                        Model const &model, bool failfast, int num_threads);

    /** Function to process warnings */
    outputFunctionType warning = printWarnMsgIdAndTxt;

    /** Function to process errors */
    outputFunctionType error = printErrMsgIdAndTxt;

    /**
     * @brief printf interface to warning()
     * @param identifier warning identifier
     * @param format string with warning message printf-style format
     * @param ... arguments to be formatted
     */
    void warningF(const char *identifier, const char *format, ...) const;

    /**
     * @brief printf interface to error()
     * @param identifier warning identifier
     * @param format string with error message printf-style format
     * @param ... arguments to be formatted
     */
    void errorF(const char *identifier, const char *format, ...) const;

    /**
     * @brief Checks the values in an array for NaNs and Infs
     *
     * @param array array
     * @param fun name of calling function
     * @return AMICI_RECOVERABLE_ERROR if a NaN/Inf value was found,
     * AMICI_SUCCESS otherwise
     */
    int checkFinite(gsl::span<const realtype> array, const char *fun);
};

/**
 * @brief Core integration routine. Initializes the solver and runs the forward
 * and backward problem.
 *
 * @param solver Solver instance
 * @param edata pointer to experimental data object
 * @param model model specification object
 * @param rethrow rethrow integration exceptions?
 * @return rdata pointer to return data object
 */
std::unique_ptr<ReturnData> runAmiciSimulation(Solver &solver,
                                               const ExpData *edata,
                                               Model &model,
                                               bool rethrow = false);

/**
 * @brief Same as runAmiciSimulation, but for multiple ExpData instances. When
 * compiled with OpenMP support, this function runs multi-threaded.
 *
 * @param solver Solver instance
 * @param edatas experimental data objects
 * @param model model specification object
 * @param failfast flag to allow early termination
 * @param num_threads number of threads for parallel execution
 * @return vector of pointers to return data objects
 */
std::vector<std::unique_ptr<ReturnData>>
runAmiciSimulations(Solver const &solver, const std::vector<ExpData *> &edatas,
                    Model const &model, bool failfast, int num_threads);

} // namespace amici

#endif /* amici_h */
