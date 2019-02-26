#ifndef AMICI_MISC_H
#define AMICI_MISC_H

#include "amici/defines.h"
#include <sunmatrix/sunmatrix_sparse.h> // SUNMatrixContent_Sparse

#include <algorithm>
#include <vector>
#include <memory>

namespace amici {

/** Checks the values in an array for NaNs and Infs
 *
 * @param array array
 * @param fun name of calling function
 * @return AMICI_RECOVERABLE_ERROR if a NaN/Inf value was found, AMICI_SUCCESS otherwise
 */
int checkFinite(std::vector <realtype> const& array, const char* fun);

/** Checks the values in an array for NaNs and Infs
 *
 * @param N number of elements in array
 * @param array array
 * @param fun name of calling function
 * @return AMICI_RECOVERABLE_ERROR if a NaN/Inf value was found, AMICI_SUCCESS otherwise
 */
int checkFinite(const int N, const realtype *array, const char* fun);


/**
  * @brief Remove parameter scaling according to the parameter scaling in pscale
  *
  * @param bufferScaled scaled parameters
  * @param pscale parameter scaling
  * @param n number of elements in bufferScaled, pscale and bufferUnscaled
  * @param bufferUnscaled unscaled parameters are written to the array
  *
  * @return status flag indicating success of execution @type int
  */
void unscaleParameters(const double *bufferScaled,
                       const ParameterScaling *pscale,
                       int n,
                       double *bufferUnscaled);

/**
  * @brief Remove parameter scaling according to the parameter scaling in pscale.
  *
  * All vectors must be of same length
  *
  * @param bufferScaled scaled parameters
  * @param pscale parameter scaling
  * @param bufferUnscaled unscaled parameters are written to the array
  */
void unscaleParameters(std::vector<double> const& bufferScaled,
                       std::vector<ParameterScaling> const& pscale,
                       std::vector<double> & bufferUnscaled);

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
void scaleParameters(const std::vector<double> &bufferUnscaled,
                     const std::vector<ParameterScaling> &pscale,
                     std::vector<double> &bufferScaled);

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

