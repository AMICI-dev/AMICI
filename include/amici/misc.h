#ifndef AMICI_MISC_H
#define AMICI_MISC_H

#include "amici/defines.h"
#include <sunmatrix/sunmatrix_sparse.h> // SUNMatrixContent_Sparse

#include <algorithm>
#include <vector>
#include <memory>
#include <regex>

#include <gsl/gsl-lite.hpp>

namespace amici {

/**
 * @brief creates a slice from existing data
 *
 * @param data to be sliced
 * @param index slice index
 * @param size slice size
 * @return span of the slice
 */

 gsl::span<realtype> slice(std::vector<realtype> &data, int index,
                           unsigned size);



/**
  * @brief Remove parameter scaling according to the parameter scaling in pscale
  *
  * All vectors must be of same length.
  *
  * @param bufferScaled scaled parameters
  * @param pscale parameter scaling
  * @param bufferUnscaled unscaled parameters are written to the array
  */
void unscaleParameters(gsl::span<const realtype> bufferScaled,
                       gsl::span<const ParameterScaling> pscale,
                       gsl::span<realtype> bufferUnscaled);

/**
  * @brief Remove parameter scaling according to `scaling`
  *
  * @param scaledParameter scaled parameter
  * @param scaling parameter scaling
  *
  * @return Unscaled parameter
  */
double getUnscaledParameter(double scaledParameter, ParameterScaling scaling);


/**
 * @brief Apply parameter scaling according to `scaling`
 * @param unscaledParameter
 * @param scaling parameter scaling
 * @return Scaled parameter
 */
double getScaledParameter(double unscaledParameter, ParameterScaling scaling);


/**
 * @brief Apply parameter scaling according to `scaling`
 * @param bufferUnscaled
 * @param pscale parameter scaling
 * @param bufferScaled destination
 */
void scaleParameters(gsl::span<const realtype> bufferUnscaled,
                     gsl::span<const ParameterScaling> pscale,
                     gsl::span<realtype> bufferScaled);

/**
 * @brief Returns the current backtrace as std::string
 * @param maxFrames Number of frames to include
 * @return Backtrace
 */
std::string backtraceString(int maxFrames);

/**
 * @brief Convert std::regex_constants::error_type to string
 * @param err_type error type
 * @return Error type as string
 */
std::string regexErrorToString(std::regex_constants::error_type err_type);

/**
 * @brief Format printf-style arguments to std::string
 * @param fmt Format string
 * @param ap Argument list pointer
 * @return Formatted String
 */
std::string printfToString(const char *fmt, va_list ap);


} // namespace amici

#endif // AMICI_MISC_H
