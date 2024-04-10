#include "amici/misc.h"
#include "amici/symbolic_functions.h"

#include <cstdarg>
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
#include <cxxabi.h> // for __cxa_demangle
#include <dlfcn.h>  // for dladdr
#include <execinfo.h>
#endif

namespace amici {

void writeSlice(AmiVector const& s, gsl::span<realtype> b) {
    writeSlice(s.getVector(), b);
};

double getUnscaledParameter(double scaledParameter, ParameterScaling scaling) {
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

void unscaleParameters(
    gsl::span<realtype const> bufferScaled,
    gsl::span<ParameterScaling const> pscale, gsl::span<realtype> bufferUnscaled
) {
    Expects(bufferScaled.size() == pscale.size());
    Expects(bufferScaled.size() == bufferUnscaled.size());

    for (gsl::span<realtype>::index_type ip = 0; ip < bufferScaled.size();
         ++ip) {
        bufferUnscaled[ip] = getUnscaledParameter(bufferScaled[ip], pscale[ip]);
    }
}

double getScaledParameter(double unscaledParameter, ParameterScaling scaling) {
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

void scaleParameters(
    gsl::span<realtype const> bufferUnscaled,
    gsl::span<ParameterScaling const> pscale, gsl::span<realtype> bufferScaled
) {
    Expects(bufferScaled.size() == pscale.size());
    Expects(bufferScaled.size() == bufferUnscaled.size());

    for (gsl::span<realtype>::index_type ip = 0; ip < bufferUnscaled.size();
         ++ip) {
        bufferScaled[ip] = getScaledParameter(bufferUnscaled[ip], pscale[ip]);
    }
}

std::string backtraceString(int const maxFrames, int const first_frame) {
    std::ostringstream trace_buf;

#ifdef PLATFORM_WINDOWS
    trace_buf << "stacktrace not available on windows platforms\n";
#else
    int const last_frame = first_frame + maxFrames;
    void* callstack[last_frame];
    char buf[1024];
    int nFrames = backtrace(callstack, last_frame);
    char** symbols = backtrace_symbols(callstack, nFrames);

    for (int i = first_frame; i < nFrames; i++) {
        // call
        Dl_info info;
        if (dladdr(callstack[i], &info) && info.dli_sname) {
            char* demangled = nullptr;
            int status = -1;
            if (info.dli_sname[0] == '_')
                demangled = abi::__cxa_demangle(
                    info.dli_sname, nullptr, nullptr, &status
                );
            snprintf(
                buf, sizeof(buf), "%-3d %*p %s + %zd\n", i - 2,
                int(2 + sizeof(void*) * 2), callstack[i],
                status == 0                 ? demangled
                : info.dli_sname == nullptr ? symbols[i]
                                            : info.dli_sname,
                static_cast<ssize_t>(
                    (char*)callstack[i] - (char*)info.dli_saddr
                )
            );
            free(demangled);
        } else {
            snprintf(
                buf, sizeof(buf), "%-3d %*p %s\n", i - 2,
                int(2 + sizeof(void*) * 2), callstack[i], symbols[i]
            );
        }
        trace_buf << buf;
    }
    free(symbols);

    if (nFrames == last_frame)
        trace_buf << "[possibly truncated]\n";
#endif
    return trace_buf.str();
}

std::string regexErrorToString(std::regex_constants::error_type err_type) {
    switch (err_type) {
    case std::regex_constants::error_collate:
        return "error_collate";
    case std::regex_constants::error_ctype:
        return "error_ctype";
    case std::regex_constants::error_escape:
        return "error_escape";
    case std::regex_constants::error_backref:
        return "error_backref";
    case std::regex_constants::error_brack:
        return "error_brack";
    case std::regex_constants::error_paren:
        return "error_paren";
    case std::regex_constants::error_brace:
        return "error_brace";
    case std::regex_constants::error_badbrace:
        return "error_badbrace";
    case std::regex_constants::error_range:
        return "error_range";
    case std::regex_constants::error_space:
        return "error_space";
    case std::regex_constants::error_badrepeat:
        return "error_badrepeat";
    case std::regex_constants::error_complexity:
        return "error_complexity";
    case std::regex_constants::error_stack:
        return "error_stack";
    default:
        return "unknown error";
    }
}

std::string printfToString(char const* fmt, va_list ap) {
    // Get size of string
    va_list ap_count;
    va_copy(ap_count, ap);
    auto size = vsnprintf(nullptr, 0, fmt, ap_count);
    va_end(ap_count);
    ++size;

    // actual formatting
    auto buf = new char[size];
    size = vsnprintf(buf, size, fmt, ap);
    std::string str(buf, size);
    delete[] buf;

    return str;
}

std::pair<size_t, size_t> unravel_index(size_t flat_idx, size_t num_cols) {
    return {flat_idx / num_cols, flat_idx % num_cols};
}

} // namespace amici
