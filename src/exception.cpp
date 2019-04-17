#include "amici/exception.h"

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
#include <execinfo.h>
#include <dlfcn.h>    // for dladdr
#include <cxxabi.h>   // for __cxa_demangle
#endif

namespace amici {

AmiException::AmiException(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(msg, sizeof(msg), fmt, ap);
    va_end(ap);
    storeBacktrace(12);
}

AmiException::AmiException(const amici::AmiException &old) {
    snprintf(msg, sizeof(msg), "%s", old.msg);
    snprintf(trace, sizeof(trace), "%s", old.trace);
}

const char *AmiException::what() const noexcept {
    return msg;
}

const char *AmiException::getBacktrace() const {
    return trace;
}

void AmiException::storeBacktrace(const int nMaxFrames) {
    std::ostringstream trace_buf;

#ifdef PLATFORM_WINDOWS

    trace_buf << "stacktrace not available on windows platforms\n";

#else

    void *callstack[nMaxFrames];
    char buf[1024];
    int nFrames = backtrace(callstack, nMaxFrames);
    char **symbols = backtrace_symbols(callstack, nFrames);

    for (int i = 2; i < nFrames;
         i++) { // start at 2 to omit AmiException and storeBacktrace
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
    if (nFrames == nMaxFrames)
        trace_buf << "[truncated]\n";

#endif

    snprintf(trace, sizeof(trace), "%s", trace_buf.str().c_str());
}

CvodeException::CvodeException(const int error_code, const char *function) :
    AmiException("Cvode routine %s failed with error code %i",function,error_code){}

IDAException::IDAException(const int error_code, const char *function) :
    AmiException("IDA routine %s failed with error code %i",function,error_code){}

IntegrationFailure::IntegrationFailure(int code, realtype t) :
    AmiException("AMICI failed to integrate the forward problem"),
    error_code(code), time(t) {}

IntegrationFailureB::IntegrationFailureB(int code, realtype t) :
    AmiException("AMICI failed to integrate the backward problem"),
    error_code(code), time(t) {}

NewtonFailure::NewtonFailure(int code, const char *function) :
    AmiException("NewtonSolver routine %s failed with error code %i",function,code) {
    error_code = code;
}

} // namespace amici
