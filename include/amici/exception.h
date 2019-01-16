#ifndef amici_exception_h
#define amici_exception_h

#include "amici/defines.h" // necessary for realtype

#include <exception>
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

    /** @brief amici exception handler class
     *
     * has a printf style interface to allow easy generation of error messages
     */
    class AmiException : public std::exception {
    public:
        AmiException(char const* fmt, ...) : std::exception() {
            /** constructor with printf style interface
             * @param fmt error message with printf format
             * @param ... printf formating variables
             */
            va_list ap;
            va_start(ap, fmt);
            vsnprintf(msg, sizeof(msg), fmt, ap);
            va_end(ap);
            storeBacktrace(12);
        }

        AmiException(const AmiException& old) {
            /** copy constructor
             * @param old object to copy from
             */
            snprintf(msg, sizeof(msg), "%s", old.msg);
            snprintf(trace, sizeof(trace), "%s", old.trace);
        }

        /** override of default error message function
         * @return msg error message
         */
        const char* what() const throw() {
            return msg;
        }

        /** returns the stored backtrace
         * @return trace backtrace
         */
        const char *getBacktrace() const {
            return trace;
        }

        /** stores the current backtrace
         * @param nMaxFrames number of frams to go back in stacktrace
         */
        void storeBacktrace(const int nMaxFrames) {
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
                    char *demangled = NULL;
                    int status = -1;
                    if (info.dli_sname[0] == '_')
                        demangled = abi::__cxa_demangle(info.dli_sname, NULL, 0,
                                                        &status);
                    snprintf(buf, sizeof(buf), "%-3d %*p %s + %zd\n", i - 2,
                             int(2 + sizeof(void *) * 2), callstack[i],
                             status == 0 ? demangled
                                         : info.dli_sname == 0 ? symbols[i]
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

    private:
        char msg[500];
        char trace[500];
    };

    /** @brief cvode exception handler class
     */
    class CvodeException : public AmiException  {
    public:
        /** constructor
         * @param error_code error code returned by cvode function
         * @param function cvode function name
         */
        CvodeException(const int error_code, const char *function) :
        AmiException("Cvode routine %s failed with error code %i",function,error_code){}
    };

    /** @brief ida exception handler class
     */
    class IDAException : public AmiException  {
    public:
        /** constructor
         * @param error_code error code returned by ida function
         * @param function ida function name
         */
        IDAException(const int error_code, const char *function) :
        AmiException("IDA routine %s failed with error code %i",function,error_code){}
    };

    /** @brief integration failure exception for the forward problem
     * this exception should be thrown when an integration failure occured
     * for this exception we can assume that we can recover from the exception
     * and return a solution struct to the user
     */
    class IntegrationFailure : public AmiException  {
    public:
        /** error code returned by cvode/ida */
        int error_code;
        /** time of integration failure */
        realtype time;
        /** constructor
         * @param code error code returned by cvode/ida
         * @param t time of integration failure
         */
        IntegrationFailure(int code, realtype t) :
        AmiException("AMICI failed to integrate the forward problem"),
        error_code(code), time(t) {}
    };

    /** @brief integration failure exception for the backward problem
     * this exception should be thrown when an integration failure occured
     * for this exception we can assume that we can recover from the exception
     * and return a solution struct to the user
     */
    class IntegrationFailureB : public AmiException  {
    public:
        /** error code returned by cvode/ida */
        int error_code;
        /** time of integration failure */
        realtype time;
        /** constructor
         * @param code error code returned by cvode/ida
         * @param t time of integration failure
         */
        IntegrationFailureB(int code, realtype t) :
        AmiException("AMICI failed to integrate the backward problem"),
        error_code(code), time(t) {}
    };

    /** @brief setup failure exception
     * this exception should be thrown when the solver setup failed
     * for this exception we can assume that we cannot recover from the exception
     * and an error will be thrown
     */
    class SetupFailure : public AmiException {
    public:
        /** constructor, simply calls AmiException constructor
         * @param msg
         */
        explicit SetupFailure(const char *msg) : AmiException(msg) {}
    };

    /** @brief newton failure exception
     * this exception should be thrown when the steady state computation
     * failed to converge for this exception we can assume that we can
     * recover from the exception and return a solution struct to the user
     */
    class NewtonFailure : public AmiException {
    public:
        /** error code returned by solver */
        int error_code;
        /** constructor, simply calls AmiException constructor
         * @param function name of the function in which the error occured
         * @param code error code
         */
        NewtonFailure(int code, const char *function) :
        AmiException("NewtonSolver routine %s failed with error code %i",function,code) {
            error_code = code;
        }
    };


}

#endif /* amici_exception_h */
