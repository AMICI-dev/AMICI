#include <amici/defines.h>
#include <array>

namespace amici {

namespace model_model_robertson_py {

// clang-format off

std::array<const char*, 3> parameterNames = {
    "p1", // p[0]
"p2", // p[1]
"p3", // p[2]
};

std::array<const char*, 1> fixedParameterNames = {
    "k1", // k[0]
};

std::array<const char*, 3> stateNames = {
    "x1", // x_rdata[0]
"x2", // x_rdata[1]
"x3", // x_rdata[2]
};

std::array<const char*, 3> observableNames = {
    "y0", // y[0]
"y1", // y[1]
"y2", // y[2]
};

std::array<const ObservableScaling, 3> observableScalings = {
    ObservableScaling::lin, // y[0]
ObservableScaling::lin, // y[1]
ObservableScaling::lin, // y[2]
};

std::array<const char*, 1> expressionNames = {
    "flux_r0", // w[0]
};

std::array<const char*, 3> parameterIds = {
    "p1", // p[0]
"p2", // p[1]
"p3", // p[2]
};

std::array<const char*, 1> fixedParameterIds = {
    "k1", // k[0]
};

std::array<const char*, 3> stateIds = {
    "x1", // x_rdata[0]
"x2", // x_rdata[1]
"x3", // x_rdata[2]
};

std::array<const char*, 3> observableIds = {
    "obs_x1", // y[0]
"obs_x2", // y[1]
"obs_x3", // y[2]
};

std::array<const char*, 1> expressionIds = {
    "flux_r0", // w[0]
};

std::array<int, 3> stateIdxsSolver = {
    0, 1, 2
};

// clang-format on

} // namespace model_model_robertson_py

} // namespace amici
