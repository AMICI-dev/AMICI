#include "amici/misc.h"
#include "amici/symbolic_functions.h"

#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <ctime>
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

#if defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h> // for GetThreadTimes
#else
#include <time.h> // for clock_gettime
#endif

namespace amici {

void write_slice(AmiVector const& s, gsl::span<realtype> const b) {
    write_slice(s.get_vector(), b);
}

double unscale_parameter(
    double const scaledParameter, ParameterScaling const scaling
) {
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

void unscale_parameters(
    gsl::span<realtype const> const bufferScaled,
    gsl::span<ParameterScaling const> const pscale,
    gsl::span<realtype> const bufferUnscaled
) {
    Expects(bufferScaled.size() == pscale.size());
    Expects(bufferScaled.size() == bufferUnscaled.size());

    for (gsl::span<realtype>::index_type ip = 0; ip < bufferScaled.size();
         ++ip) {
        bufferUnscaled[ip] = unscale_parameter(bufferScaled[ip], pscale[ip]);
    }
}

double scale_parameter(
    double const unscaledParameter, ParameterScaling const scaling
) {
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

void scale_parameters(
    gsl::span<realtype const> const bufferUnscaled,
    gsl::span<ParameterScaling const> const pscale,
    gsl::span<realtype> const bufferScaled
) {
    Expects(bufferScaled.size() == pscale.size());
    Expects(bufferScaled.size() == bufferUnscaled.size());

    for (gsl::span<realtype>::index_type ip = 0; ip < bufferUnscaled.size();
         ++ip) {
        bufferScaled[ip] = scale_parameter(bufferUnscaled[ip], pscale[ip]);
    }
}

std::string get_backtrace_string(int const maxFrames, int const first_frame) {
    std::ostringstream trace_buf;

#ifdef PLATFORM_WINDOWS
    trace_buf << "stacktrace not available on windows platforms\n";
#else
    int const last_frame = first_frame + maxFrames;
    std::vector<void*> callstack(last_frame);
    char buf[1024];
    int const nFrames = backtrace(callstack.data(), last_frame);
    char** symbols = backtrace_symbols(callstack.data(), nFrames);

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

std::string
regex_error_to_string(std::regex_constants::error_type const err_type) {
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

std::string printf_to_string(char const* fmt, va_list ap) {
    // Get size of string
    va_list ap_count;
    va_copy(ap_count, ap);
    auto size = vsnprintf(nullptr, 0, fmt, ap_count);
    va_end(ap_count);
    ++size;

    // actual formatting
    auto const buf = new char[size];
    size = vsnprintf(buf, size, fmt, ap);
    std::string str(buf, size);
    delete[] buf;

    return str;
}

std::pair<size_t, size_t>
unravel_index(size_t const flat_idx, size_t const num_cols) {
    return {flat_idx / num_cols, flat_idx % num_cols};
}

bool const CpuTimer::uses_thread_clock =
#if defined(_WIN32) || defined(CLOCK_THREAD_CPUTIME_ID)
    true;
#else
    false;
#endif

double CpuTimer::now_seconds() {
#if defined(_WIN32)
    FILETIME creation_time, exit_time, kernel_time, user_time;
    GetThreadTimes(
        GetCurrentThread(), &creation_time, &exit_time, &kernel_time, &user_time
    );
    // FILETIME counts 100-nanosecond intervals
    auto to_seconds = [](FILETIME const& ft) {
        auto const ticks = (static_cast<std::uint64_t>(ft.dwHighDateTime) << 32)
                           | ft.dwLowDateTime;
        return static_cast<double>(ticks) * 1e-7;
    };
    return to_seconds(kernel_time) + to_seconds(user_time);
#elif defined(CLOCK_THREAD_CPUTIME_ID)
    timespec ts{};
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts);
    return static_cast<double>(ts.tv_sec)
           + static_cast<double>(ts.tv_nsec) * 1e-9;
#else
    return static_cast<double>(std::clock()) / CLOCKS_PER_SEC;
#endif
}

} // namespace amici
