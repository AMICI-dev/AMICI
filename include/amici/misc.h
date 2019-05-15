#ifndef AMICI_MISC_H
#define AMICI_MISC_H

#include "amici/defines.h"
#include <sunmatrix/sunmatrix_sparse.h> // SUNMatrixContent_Sparse

#include <algorithm>
#include <vector>
#include <memory>

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
 
 gsl::span<realtype> slice(std::vector<realtype> &data, const int index,
                           const unsigned size);

/**
 * @brief Checks the values in an array for NaNs and Infs
 *
 * @param array array
 * @param fun name of calling function
 * @return AMICI_RECOVERABLE_ERROR if a NaN/Inf value was found, AMICI_SUCCESS otherwise
 */
int checkFinite(gsl::span<const realtype> array, const char* fun);


/**
  * @brief Remove parameter scaling according to the parameter scaling in pscale
  *
  * All vectors must be of same length.
  *
  * @param bufferScaled scaled parameters
  * @param pscale parameter scaling
  * @param bufferUnscaled unscaled parameters are written to the array
  *
  * @return status flag indicating success of execution @type int
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


} // namespace amici

#ifndef __cpp_lib_make_unique
// custom make_unique while we are still using c++11
namespace std {
template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}
}
#endif

#endif // AMICI_MISC_H

