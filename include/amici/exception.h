#ifndef amici_exception_h
#define amici_exception_h

#include "amici/defines.h" // necessary for realtype

#include <exception>
#include <array>

namespace amici {

/**
 * @brief AMICI exception class
 *
 * Has a printf style interface to allow easy generation of error messages
 */
class AmiException : public std::exception {
public:
    /**
     * @brief Constructor with printf style interface
     * @param fmt error message with printf format
     * @param ... printf formatting variables
     */
    explicit AmiException(char const* fmt, ...);

    /**
     * @brief Override of default error message function
     * @return msg error message
     */
    const char* what() const noexcept override;

    /**
     * @brief Returns the stored backtrace
     * @return trace backtrace
     */
    const char *getBacktrace() const;

    /**
     * @brief Stores the current backtrace
     * @param nMaxFrames number of frames to go back in stacktrace
     */
    void storeBacktrace(int nMaxFrames);

private:
    std::array<char, 500> msg_;
    std::array<char, 500> trace_;
};


/**
 * @brief cvode exception handler class
 */
class CvodeException : public AmiException  {
public:
    /**
     * @brief Constructor
     * @param error_code error code returned by cvode function
     * @param function cvode function name
     */
    CvodeException(int error_code, const char *function);
};


/**
 * @brief ida exception handler class
 */
class IDAException : public AmiException  {
public:
    /**
     * @brief Constructor
     * @param error_code error code returned by ida function
     * @param function ida function name
     */
    IDAException(int error_code, const char *function);
};


/**
 * @brief Integration failure exception for the forward problem
 *
 * This exception should be thrown when an integration failure occurred
 * for this exception we can assume that we can recover from the exception
 * and return a solution struct to the user
 */
class IntegrationFailure : public AmiException  {
  public:
    /**
     * @brief Constructor
     * @param code error code returned by cvode/ida
     * @param t time of integration failure
     */
    IntegrationFailure(int code, realtype t);

    /** error code returned by cvodes/idas */
    int error_code;

    /** time of integration failure */
    realtype time;
};


/**
 * @brief Integration failure exception for the backward problem
 *
 * This exception should be thrown when an integration failure occurred
 * for this exception we can assume that we can recover from the exception
 * and return a solution struct to the user
 */
class IntegrationFailureB : public AmiException  {
  public:
    /**
     * @brief Constructor
     * @param code error code returned by cvode/ida
     * @param t time of integration failure
     */
    IntegrationFailureB(int code, realtype t);

    /** error code returned by cvode/ida */
    int error_code;

    /** time of integration failure */
    realtype time;
};


/**
 * @brief Setup failure exception
 *
 * This exception should be thrown when the solver setup failed
 * for this exception we can assume that we cannot recover from the exception
 * and an error will be thrown
 */
class SetupFailure : public AmiException {
  public:
    using AmiException::AmiException;
};


/**
 * @brief Newton failure exception
 *
 * This exception should be thrown when the steady state computation
 * failed to converge for this exception we can assume that we can
 * recover from the exception and return a solution struct to the user
 */
class NewtonFailure : public AmiException {
public:
    /**
     * @brief Constructor, simply calls AmiException constructor
     * @param function name of the function in which the error occurred
     * @param code error code
     */
    NewtonFailure(int code, const char *function);

    /** error code returned by solver */
    int error_code;
};

} // namespace amici

#endif /* amici_exception_h */
