#include <amici/defines.h>
#include <array>

namespace amici::model_model_nested_events_py {

// clang-format off

std::array<const char*, 5> free_parameter_names = {
    "V_0", // p[0]
"V_0_inject", // p[1]
"t_0", // p[2]
"rho_V", // p[3]
"delta_V", // p[4]
};

std::array<const char*, 0> fixed_parameter_names = {
    
};

std::array<const char*, 1> state_names = {
    "Virus", // x_rdata[0]
};

std::array<const char*, 1> observable_names = {
    "y0", // y[0]
};

std::array<const ObservableScaling, 1> observable_scalings = {
    ObservableScaling::lin, // y[0]
};

std::array<const char*, 0> expression_names = {
    
};

std::array<const char*, 5> free_parameter_ids = {
    "V_0", // p[0]
"V_0_inject", // p[1]
"t_0", // p[2]
"rho_V", // p[3]
"delta_V", // p[4]
};

std::array<const char*, 0> fixed_parameter_ids = {
    
};

std::array<const char*, 1> state_ids = {
    "Virus", // x_rdata[0]
};

std::array<const char*, 1> observable_ids = {
    "obs_Virus", // y[0]
};

std::array<const char*, 0> expression_ids = {
    
};

std::array<int, 1> state_idxs_solver = {
    0
};

// clang-format on

} // namespace amici::model_model_nested_events_py
