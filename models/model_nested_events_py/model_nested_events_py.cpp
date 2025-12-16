#include <amici/defines.h>
#include <array>
#include <string_view>

namespace amici::model_model_nested_events_py {

// clang-format off

extern const std::array<std::string_view const, 5> free_parameter_names = {
    "V_0", // p[0]
"V_0_inject", // p[1]
"t_0", // p[2]
"rho_V", // p[3]
"delta_V", // p[4]
};

extern const std::array<std::string_view const, 0> fixed_parameter_names = {
    
};

extern const std::array<std::string_view const, 1> state_names = {
    "Virus", // x_rdata[0]
};

extern const std::array<std::string_view const, 1> state_names_solver = {
    "Virus", // x_solver[0]
};

extern const std::array<std::string_view const, 1> observable_names = {
    "y0", // y[0]
};

std::array<const ObservableScaling, 1> observable_scalings = {
    ObservableScaling::lin, // y[0]
};

extern const std::array<std::string_view const, 0> expression_names = {
    
};

extern const std::array<std::string_view const, 5> free_parameter_ids = {
    "V_0", // p[0]
"V_0_inject", // p[1]
"t_0", // p[2]
"rho_V", // p[3]
"delta_V", // p[4]
};

extern const std::array<std::string_view const, 0> fixed_parameter_ids = {
    
};

extern const std::array<std::string_view const, 1> state_ids = {
    "Virus", // x_rdata[0]
};

extern const std::array<std::string_view const, 1> state_ids_solver = {
    "Virus", // x_solver[0]
};

extern const std::array<std::string_view const, 1> observable_ids = {
    "obs_Virus", // y[0]
};

extern const std::array<std::string_view const, 0> expression_ids = {
    
};

std::array<int, 1> state_idxs_solver = {
    0
};

// clang-format on

} // namespace amici::model_model_nested_events_py
