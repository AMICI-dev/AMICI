#ifndef amici_exception_h
#define amici_exception_h

#include <amici_defines.h> // necessary for realtype
#include <exception>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <sstream>
#define PLATFORM_WINDOWS 1
#define PLATFORM_UNIX 2
#if defined(_WIN32)
    #define PLATFORM PLATFORM_WINDOWS // Windows
#elif defined(_WIN64)
    #define PLATFORM PLATFORM_WINDOWS // Windows
#elif defined(__CYGWIN__) && !defined(_WIN32)
    #define PLATFORM PLATFORM_WINDOWS // Windows (Cygwin POSIX under Microsoft Window)
#else
    #define PLATFORM PLATFORM_UNIX
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
             * @param[in] fmt error message with printf format
             * @param[in] ... printf formating variables
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
            
            #if PLATFORM == PLATFORM_WINDOWS
            
            trace_buf << "stacktrace not available on windows platforms\n";
            
            #else
            
            void *callstack[nMaxFrames];
            char buf[1024];
            int nFrames = backtrace(callstack, nMaxFrames);
            char **symbols = backtrace_symbols(callstack, nFrames);
            
            for (int i = 0; i < nFrames; i++) {
                Dl_info info;
                if (dladdr(callstack[i], &info) && info.dli_sname) {
                    char *demangled = NULL;
                    int status = -1;
                    if (info.dli_sname[0] == '_')
                        demangled = abi::__cxa_demangle(info.dli_sname, NULL, 0, &status);
                    snprintf(buf, sizeof(buf), "%-3d %*p %s + %zd\n",
                             i, int(2 + sizeof(void*) * 2), callstack[i],
                             status == 0 ? demangled :
                             info.dli_sname == 0 ? symbols[i] : info.dli_sname,
                             (char *)callstack[i] - (char *)info.dli_saddr);
                    free(demangled);
                } else {
                    snprintf(buf, sizeof(buf), "%-3d %*p %s\n",
                             i, int(2 + sizeof(void*) * 2), callstack[i], symbols[i]);
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
         * @param[in] error_code error code returned by cvode function
         * @param[in] function cvode function name
         */
        CvodeException(const int error_code, const char *function) :
        AmiException("Cvode routine %s failed with error code (%i)",function,error_code){}
    };
    
    /** @brief ida exception handler class
     */
    class IDAException : public AmiException  {
    public:
        /** constructor
         * @param[in] error_code error code returned by ida function
         * @param[in] function ida function name
         */
        IDAException(const int error_code, const char *function) :
        AmiException("IDA routine %s failed with error code (%i)",function,error_code){}
    };
    
    /** @brief integration failure exception
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
         * @param[in] code error code returned by cvode/ida
         * @param[in] t time of integration failure
         */
        IntegrationFailure(int code, realtype t) :
        AmiException("AMICI failed to integrate the problem") {
            error_code = code;
            time = t;
        }
    };
    
    /** @brief setup failure exception
     * this exception should be thrown when the solver setup failed
     * for this exception we can assume that we cannot recover from the exception
     * and an error will be thrown
     */
    class SetupFailure : public AmiException {
    public:
        /** constructor, simply calls AmiException constructor
         * @param[in] msg
         */
        SetupFailure(const char *msg) : AmiException(msg) {}
    };
    
    /** @brief newton failure exception
     * this exception should be thrown when the steady state computation
     * failed to converge for this exception we can assume that we can
     * recover from the exception and return a solution struct to the user
     */
    class NewtonFailure : public AmiException {
    public:
        /** constructor, simply calls AmiException constructor
         * @param[in] msg
         */
        NewtonFailure(const char *msg) : AmiException(msg) {}
    };
    
    
}

#endif /* amici_exception_h */
