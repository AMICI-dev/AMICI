#include <amici/defines.h>
#include <array>

namespace amici {

namespace model_model_nested_events_py {

// clang-format off

std::array<const char*, 5> parameterNames = {
    "V_0", // p[0]
"V_0_inject", // p[1]
"t_0", // p[2]
"rho_V", // p[3]
"delta_V", // p[4]
};

std::array<const char*, 0> fixedParameterNames = {
    
};

std::array<const char*, 1> stateNames = {
    "Virus", // x_rdata[0]
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

std::array<const char*, 5> parameterIds = {
    "V_0", // p[0]
"V_0_inject", // p[1]
"t_0", // p[2]
"rho_V", // p[3]
"delta_V", // p[4]
};

std::array<const char*, 0> fixedParameterIds = {
    
};

std::array<const char*, 1> stateIds = {
    "Virus", // x_rdata[0]
};

std::array<const char*, 1> observableIds = {
    "obs_Virus", // y[0]
};

std::array<const char*, 1> expressionIds = {
    "flux_r0", // w[0]
};

std::array<int, 1> stateIdxsSolver = {
    0
};

// clang-format on

} // namespace model_model_nested_events_py

} // namespace amici
