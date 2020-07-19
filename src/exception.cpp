#include "amici/exception.h"
#include "amici/misc.h"

#include <cstdarg>
#include <cstdio>
#include <cstring>


namespace amici {

AmiException::AmiException(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(msg_.data(), msg_.size(), fmt, ap);
    va_end(ap);
    storeBacktrace(12);
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
