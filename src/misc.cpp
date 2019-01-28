#include "amici/misc.h"
#include "amici/amici.h"
#include "amici/symbolic_functions.h"

#include <cstdio>
#include <cstring>

namespace amici {

int checkFinite(const int N, const realtype *array, const char *fun) {
    for (int idx = 0; idx < N; idx++) {
        if (isNaN(array[idx])) {
            warnMsgIdAndTxt(
                "AMICI:NaN",
                "AMICI encountered a NaN value at index %i of %i in %s!", idx,
                N, fun);
            return AMICI_RECOVERABLE_ERROR;
        }
        if (isInf(array[idx])) {
            warnMsgIdAndTxt(
                "AMICI:Inf",
                "AMICI encountered an Inf value at index %i of %i in %s!", idx,
                N, fun);
            return AMICI_RECOVERABLE_ERROR;
        }
    }
    return AMICI_SUCCESS;
}

double getUnscaledParameter(double scaledParameter, ParameterScaling scaling)
{
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

void unscaleParameters(const double *bufferScaled, const ParameterScaling *pscale, int n, double *bufferUnscaled)
{
    for (int ip = 0; ip < n; ++ip) {
        bufferUnscaled[ip] = getUnscaledParameter(bufferScaled[ip], pscale[ip]);
    }

}

void unscaleParameters(const std::vector<double> &bufferScaled, const std::vector<ParameterScaling> &pscale, std::vector<double> &bufferUnscaled)
{
    if(bufferScaled.size() != pscale.size() || pscale.size() != bufferUnscaled.size())
        throw AmiException("Vector size mismatch in unscaleParameters.");
    unscaleParameters(bufferScaled.data(), pscale.data(), bufferScaled.size(), bufferUnscaled.data());
}

double getScaledParameter(double unscaledParameter, ParameterScaling scaling)
{
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


void scaleParameters(const std::vector<double> &bufferUnscaled, const std::vector<ParameterScaling> &pscale, std::vector<double> &bufferScaled)
{
    if(bufferScaled.size() != pscale.size() || pscale.size() != bufferUnscaled.size())
        throw AmiException("Vector size mismatch in scaleParameters.");
    for (int ip = 0; ip < (int) bufferUnscaled.size(); ++ip) {
        bufferScaled[ip] = getScaledParameter(bufferUnscaled[ip], pscale[ip]);
    }

}

int checkFinite(const std::vector<realtype> &array, const char *fun)
{
    return checkFinite(array.size(), array.data(), fun);
}

} // namespace amici
