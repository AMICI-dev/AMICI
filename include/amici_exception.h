#ifndef amici_exception_h
#define amici_exception_h

#include <exception>
#include <cstdarg>
#include <cstdio>

namespace amici {
    
#define AMICI_ERROR_UDATA              -99
#define AMICI_ERROR_EDATA              -98
#define AMICI_ERROR_RDATA              -97
#define AMICI_ERROR_TDATA              -96
#define AMICI_ERROR_SETUP              -95
#define AMICI_ERROR_SETUPB             -94
#define AMICI_ERROR_NOTHINGTODO        -93
#define AMICI_ERROR_FSA                -92
#define AMICI_ERROR_ASA                -91
#define AMICI_ERROR_SA                 -90
#define AMICI_ERROR_SS_SENSIS          -89
#define AMICI_ERROR_DATA               -88
#define AMICI_ERROR_EVENT              -87
#define AMICI_ERROR_SIMULATION         -86
#define AMICI_ERROR_NEWTONSOLVER       -85
#define AMICI_ERROR_NEWTON_LINSOLVER   -84
#define AMICI_ERROR_NOT_IMPLEMENTED    -83
#define AMICI_ERROR_MODEL              -82
#define AMICI_ERROR_OTHER              -81
#define AMICI_ERROR_SIM2STEADYSTATE    -80
#define AMICI_ERROR_PREEQUILIBRATION   -79
#define AMICI_SUCCESS                    0

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
