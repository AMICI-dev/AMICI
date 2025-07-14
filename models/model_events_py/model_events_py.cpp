#include <amici/defines.h>
#include <array>

namespace amici {

namespace model_model_events_py {

// clang-format off

std::array<const char*, 4> parameterNames = {
    "p1", // p[0]
"p2", // p[1]
"p3", // p[2]
"p4", // p[3]
};

std::array<const char*, 4> fixedParameterNames = {
    "k1", // k[0]
"k2", // k[1]
"k3", // k[2]
"k4", // k[3]
};

std::array<const char*, 3> stateNames = {
    "x1", // x_rdata[0]
"x2", // x_rdata[1]
"x3", // x_rdata[2]
};

std::array<const char*, 1> observableNames = {
    "y1", // y[0]
};

std::array<const ObservableScaling, 1> observableScalings = {
    ObservableScaling::lin, // y[0]
};

std::array<const char*, 1> expressionNames = {
    "flux_r0", // w[0]
};

std::array<const char*, 4> parameterIds = {
    "p1", // p[0]
"p2", // p[1]
"p3", // p[2]
"p4", // p[3]
};

std::array<const char*, 4> fixedParameterIds = {
    "k1", // k[0]
"k2", // k[1]
"k3", // k[2]
"k4", // k[3]
};

std::array<const char*, 3> stateIds = {
    "x1", // x_rdata[0]
"x2", // x_rdata[1]
"x3", // x_rdata[2]
};

std::array<const char*, 1> observableIds = {
    "y1", // y[0]
};

std::array<const char*, 1> expressionIds = {
    "flux_r0", // w[0]
};

std::array<int, 3> stateIdxsSolver = {
    0, 1, 2
};

// clang-format on

} // namespace model_model_events_py

} // namespace amici
