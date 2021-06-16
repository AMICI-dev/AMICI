#ifndef AMICI_MISC_H
#define AMICI_MISC_H

#include "amici/defines.h"
#include "amici/exception.h"
#include "amici/vector.h"
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

template <class T>
gsl::span<T> slice(std::vector<T> &data, int index, unsigned size) {
    if ((index + 1) * size > data.size())
        throw std::out_of_range("requested slice is out of data range");
    if (size > 0)
        return gsl::make_span(&data.at(index*size), size);

    return gsl::make_span(static_cast<T*>(nullptr), 0);
}

/**
 * @brief creates a constant slice from existing constant data
 *
 * @param data to be sliced
 * @param index slice index
 * @param size slice size
 * @return span of the slice
 */

template <class T>
gsl::span<const T> slice(const std::vector<T> &data,
                         int index, unsigned size) {
    if ((index + 1) * size > data.size())
        throw std::out_of_range("requested slice is out of data range");
    if (size > 0)
        return gsl::make_span(&data.at(index*size), size);

    return gsl::make_span(static_cast<T*>(nullptr), 0);
}

/**
 * @brief local helper to check whether the provided buffer has the expected
 * size
 * @param buffer buffer to which values are to be written
 * @param expected_size expected size of the buffer
 */
template <class T>
void checkBufferSize(gsl::span<T> buffer,
                     typename gsl::span<T>::index_type expected_size) {
    if (buffer.size() != expected_size)
        throw AmiException("Incorrect buffer size! Was %u, expected %u.",
                           buffer.size(), expected_size);
}

/* TODO: templating writeSlice breaks implicit conversion between vector & span
 not sure whether this is fixable */

/**
 * @brief local helper function to write computed slice to provided buffer (span)
 * @param slice computed value
 * @param buffer buffer to which values are to be written
 */
template <class T>
void writeSlice(const gsl::span<const T> slice, gsl::span<T> buffer) {
    checkBufferSize(buffer, slice.size());
    std::copy(slice.begin(), slice.end(), buffer.data());
};

/**
 * @brief local helper function to write computed slice to provided buffer (vector)
 * @param s computed value
 * @param b buffer to which values are to be written
 */
template <class T>
void writeSlice(const std::vector<T> &s, std::vector<T> &b) {
    writeSlice(gsl::make_span(s.data(), s.size()),
               gsl::make_span(b.data(), b.size()));
};

/**
 * @brief local helper function to write computed slice to provided buffer (vector/span)
 * @param s computed value
 * @param b buffer to which values are to be written
 */
template <class T>
void writeSlice(const std::vector<T> &s, gsl::span<T> b) {
    writeSlice(gsl::make_span(s.data(), s.size()), b);
};

/**
 * @brief local helper function to write computed slice to provided buffer (AmiVector/span)
 * @param s computed value
 * @param b buffer to which values are to be written
 */
void writeSlice(const AmiVector &s, gsl::span<realtype> b);


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

/**
 * @brief Generic implementation for a context manager, explicitly deletes copy
 * and move operators for derived classes
 */
class ContextManager{
  public:
    ContextManager() = default;
    ContextManager(ContextManager &other) = delete;
    ContextManager(ContextManager &&other) = delete;
};

} // namespace amici

#endif // AMICI_MISC_H
