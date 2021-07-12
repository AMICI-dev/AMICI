#include "amici/exception.h"
#include "amici/misc.h"

#include <cstdarg>
#include <cstdio>
#include <cstring>


namespace amici {

AmiException::AmiException()
{
    storeBacktrace(12);
}

AmiException::AmiException(const char *fmt, ...)
    : AmiException()
{
    va_list ap;
    va_start(ap, fmt);
    storeMessage(fmt, ap);
    va_end(ap);
}

const char *AmiException::what() const noexcept {
    return msg_.data();
}

const char *AmiException::getBacktrace() const {
    return trace_.data();
}

void AmiException::storeBacktrace(const int nMaxFrames) {
    snprintf(trace_.data(), trace_.size(), "%s",
             backtraceString(nMaxFrames).c_str());
}

void AmiException::storeMessage(const char *fmt, va_list argptr)
{
    vsnprintf(msg_.data(), msg_.size(), fmt, argptr);
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

SetupFailure::SetupFailure(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    storeMessage(fmt, ap);
    va_end(ap);
}

} // namespace amici
