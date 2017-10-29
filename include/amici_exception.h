#ifndef amici_exception_h
#define amici_exception_h

#include <exception>
#include <cstdarg>
#include <cstdio>

namespace amici {

    class AmiException : public std::exception {
    public:
        AmiException(char const* fmt, ...) __attribute__((format(printf,2,3))) {
            va_list ap;
            va_start(ap, fmt);
            vsnprintf(msg, sizeof msg, fmt, ap);
            va_end(ap);
        }
        const char* what() const throw() {
            return msg;
        }
    private:
        char msg[1000];
    };
    
    class CvodeException : public AmiException  {
    public:
        CvodeException(const int error_code, const char *function) :
        AmiException("Cvode routine %s failed with error code (%i)",function,error_code){}
    };
    
    class IDAException : public AmiException  {
    public:
        IDAException(const int error_code, const char *function) :
        AmiException("IDA routine %s failed with error code (%i)",function,error_code){}
    };
    
    class NullPointerException : public AmiException  {
    public:
        NullPointerException(const char *variable) :
        AmiException("AMICI encountered a null pointer for variable (%s)",variable){}
    };
    
    
}

#endif /* amici_exception_h */
