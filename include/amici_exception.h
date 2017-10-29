#ifndef amici_exception_h
#define amici_exception_h

#include <exception>
#include <cstdarg>
#include <cstdio>

namespace amici {

    /** @brief amici exception handler class
     *
     * has a printf style interface to allow easy generation of error messages
     */
    class AmiException : public std::exception {
    public:
        AmiException(char const* fmt, ...) {
            /** constructor with printf style interface
             * @param[in] fmt error message with printf format
             * @param[in] ... printf formating variables
             */
            va_list ap;
            va_start(ap, fmt);
            vsnprintf(msg, sizeof msg, fmt, ap);
            va_end(ap);
        }
        /** override of default error message function
         * @return msg error message
         */
        const char* what() const throw() {
            return msg;
        }
    private:
        char msg[1000];
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
    
    /** @brief null pointer exception handler class
     */
    class NullPointerException : public AmiException  {
    public:
        /** constructor
         * @param[in] variable name of variable that was supposed to be accesssed
         */
        NullPointerException(const char *variable) :
        AmiException("AMICI encountered a null pointer for variable (%s)",variable){}
    };
    
    
}

#endif /* amici_exception_h */
