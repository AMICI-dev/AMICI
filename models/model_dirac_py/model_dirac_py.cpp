#include <amici/defines.h>
#include <array>

namespace amici::model_model_dirac_py {

// clang-format off

std::array<const char*, 4> parameterNames = {
    "p1", // p[0]
"p2", // p[1]
"p3", // p[2]
"p4", // p[3]
};

std::array<const char*, 0> fixedParameterNames = {
    
};

std::array<const char*, 2> stateNames = {
    "x1", // x_rdata[0]
"x2", // x_rdata[1]
};

std::array<const char*, 1> observableNames = {
    "y0", // y[0]
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

std::array<const char*, 0> fixedParameterIds = {
    
};

std::array<const char*, 2> stateIds = {
    "x1", // x_rdata[0]
"x2", // x_rdata[1]
};

std::array<const char*, 1> observableIds = {
    "obs_x2", // y[0]
};

std::array<const char*, 1> expressionIds = {
    "flux_r0", // w[0]
};

std::array<int, 2> stateIdxsSolver = {
    0, 1
};

// clang-format on

} // namespace amici::model_model_dirac_py
