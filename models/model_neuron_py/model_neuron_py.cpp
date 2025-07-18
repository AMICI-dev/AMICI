#include <amici/defines.h>
#include <array>

namespace amici {

namespace model_model_neuron_py {

// clang-format off

std::array<const char*, 4> parameterNames = {
    "a", // p[0]
"b", // p[1]
"c", // p[2]
"d", // p[3]
};

std::array<const char*, 2> fixedParameterNames = {
    "v0", // k[0]
"I0", // k[1]
};

std::array<const char*, 2> stateNames = {
    "v", // x_rdata[0]
"u", // x_rdata[1]
};

std::array<const char*, 1> observableNames = {
    "v", // y[0]
};

std::array<const ObservableScaling, 1> observableScalings = {
    ObservableScaling::lin, // y[0]
};

std::array<const char*, 1> expressionNames = {
    "flux_r0", // w[0]
};

std::array<const char*, 4> parameterIds = {
    "a", // p[0]
"b", // p[1]
"c", // p[2]
"d", // p[3]
};

std::array<const char*, 2> fixedParameterIds = {
    "v0", // k[0]
"I0", // k[1]
};

std::array<const char*, 2> stateIds = {
    "v", // x_rdata[0]
"u", // x_rdata[1]
};

std::array<const char*, 1> observableIds = {
    "y1", // y[0]
};

std::array<const char*, 1> expressionIds = {
    "flux_r0", // w[0]
};

std::array<int, 2> stateIdxsSolver = {
    0, 1
};

// clang-format on

} // namespace model_model_neuron_py

} // namespace amici
