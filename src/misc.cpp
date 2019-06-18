#include "amici/misc.h"
#include "amici/amici.h"
#include "amici/symbolic_functions.h"

#include <cstdio>
#include <cstring>
#include <sstream>
#if defined(_WIN32)
#define PLATFORM_WINDOWS // Windows
#elif defined(_WIN64)
#define PLATFORM_WINDOWS // Windows
#elif defined(__CYGWIN__) && !defined(_WIN32)
#define PLATFORM_WINDOWS // Windows (Cygwin POSIX under Microsoft Window)
#else
#include <execinfo.h>
#include <dlfcn.h>    // for dladdr
#include <cxxabi.h>   // for __cxa_demangle
#endif

namespace amici {
    
gsl::span<realtype> slice(std::vector<realtype> &data, const int index,
                          const unsigned size) {
    if ((index + 1) * size > data.size())
        throw std::out_of_range("requested slice is out of data range");
    if (size > 0)
        return gsl::make_span(&data.at(index*size), size);
    else
        return gsl::make_span(static_cast<realtype*>(nullptr), 0);
}

double getUnscaledParameter(double scaledParameter, ParameterScaling scaling)
{
    switch (scaling) {
    case ParameterScaling::log10:
        return pow(10, scaledParameter);
    case ParameterScaling::ln:
        return exp(scaledParameter);
    case ParameterScaling::none:
        return scaledParameter;
    }

    throw AmiException("Invalid value for ParameterScaling.");
}

void unscaleParameters(gsl::span<const realtype> bufferScaled,
                       gsl::span<const ParameterScaling> pscale,
                       gsl::span<realtype> bufferUnscaled)
{
    Expects(bufferScaled.size() == pscale.size());
    Expects(bufferScaled.size() == bufferUnscaled.size());

    for (gsl::span<realtype>::index_type ip = 0;
         ip < bufferScaled.size(); ++ip) {
        bufferUnscaled[ip] = getUnscaledParameter(bufferScaled[ip], pscale[ip]);
    }
}


double getScaledParameter(double unscaledParameter, ParameterScaling scaling)
{
    switch (scaling) {
    case ParameterScaling::log10:
        return log10(unscaledParameter);
    case ParameterScaling::ln:
        return log(unscaledParameter);
    case ParameterScaling::none:
        return unscaledParameter;
    }

    throw AmiException("Invalid value for ParameterScaling.");
}


void scaleParameters(gsl::span<const realtype> bufferUnscaled,
                     gsl::span<const ParameterScaling> pscale,
                     gsl::span<realtype> bufferScaled)
{
    Expects(bufferScaled.size() == pscale.size());
    Expects(bufferScaled.size() == bufferUnscaled.size());

    for (gsl::span<realtype>::index_type ip = 0;
         ip < bufferUnscaled.size(); ++ip) {
        bufferScaled[ip] = getScaledParameter(bufferUnscaled[ip], pscale[ip]);
    }

}

int checkFinite(gsl::span<const realtype> array, const char *fun)
{
    for (int idx = 0; idx < (int) array.size(); idx++) {
        if (isNaN(array[idx])) {
            warnMsgIdAndTxt(
                "AMICI:NaN",
                "AMICI encountered a NaN value at index %i of %i in %s!", idx,
                (int) array.size(), fun);
            return AMICI_RECOVERABLE_ERROR;
        }
        if (isInf(array[idx])) {
            warnMsgIdAndTxt(
                "AMICI:Inf",
                "AMICI encountered an Inf value at index %i of %i in %s!", idx,
                (int) array.size(), fun);
            return AMICI_RECOVERABLE_ERROR;
        }
    }
    return AMICI_SUCCESS;
}

std::string backtraceString(const int maxFrames)
{
    std::ostringstream trace_buf;

#ifdef PLATFORM_WINDOWS
    trace_buf << "stacktrace not available on windows platforms\n";
#else
    void *callstack[maxFrames];
    char buf[1024];
    int nFrames = backtrace(callstack, maxFrames);
    char **symbols = backtrace_symbols(callstack, nFrames);

    // start at 2 to omit AmiException and storeBacktrace
    for (int i = 2; i < nFrames; i++) {
        // call
        Dl_info info;
        if (dladdr(callstack[i], &info) && info.dli_sname) {
            char *demangled = nullptr;
            int status = -1;
            if (info.dli_sname[0] == '_')
                demangled = abi::__cxa_demangle(info.dli_sname, nullptr, nullptr,
                                                &status);
            snprintf(buf, sizeof(buf), "%-3d %*p %s + %zd\n", i - 2,
                     int(2 + sizeof(void *) * 2), callstack[i],
                     status == 0 ? demangled
                                 : info.dli_sname == nullptr ? symbols[i]
                                                               : info.dli_sname,
                     static_cast<ssize_t>((char *)callstack[i] -
                                          (char *)info.dli_saddr));
            free(demangled);
        } else {
            snprintf(buf, sizeof(buf), "%-3d %*p %s\n", i - 2,
                     int(2 + sizeof(void *) * 2), callstack[i],
                     symbols[i]);
        }
        trace_buf << buf;
    }
    free(symbols);

    if (nFrames == maxFrames)
        trace_buf << "[truncated]\n";
#endif
    return trace_buf.str();
}

} // namespace amici
