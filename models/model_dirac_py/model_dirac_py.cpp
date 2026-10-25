#include <amici/defines.h>
#include <array>
#include <string_view>

namespace amici::model_model_dirac_py {

// clang-format off

extern const std::array<std::string_view const, 4> free_parameter_names = {
    "p1", // p[0]
"p2", // p[1]
"p3", // p[2]
"p4", // p[3]
};

extern const std::array<std::string_view const, 0> fixed_parameter_names = {
    
};

extern const std::array<std::string_view const, 2> state_names = {
    "x1", // x_rdata[0]
"x2", // x_rdata[1]
};

extern const std::array<std::string_view const, 2> state_names_solver = {
    "x1", // x_solver[0]
"x2", // x_solver[1]
};

extern const std::array<std::string_view const, 1> observable_names = {
    "y0", // y[0]
};

std::array<const ObservableScaling, 1> observable_scalings = {
    ObservableScaling::lin, // y[0]
};

extern const std::array<std::string_view const, 0> expression_names = {
    
};

extern const std::array<std::string_view const, 4> free_parameter_ids = {
    "p1", // p[0]
"p2", // p[1]
"p3", // p[2]
"p4", // p[3]
};

extern const std::array<std::string_view const, 0> fixed_parameter_ids = {
    
};

extern const std::array<std::string_view const, 2> state_ids = {
    "x1", // x_rdata[0]
"x2", // x_rdata[1]
};

extern const std::array<std::string_view const, 2> state_ids_solver = {
    "x1", // x_solver[0]
"x2", // x_solver[1]
};

extern const std::array<std::string_view const, 1> observable_ids = {
    "obs_x2", // y[0]
};

extern const std::array<std::string_view const, 0> expression_ids = {
    
};

std::array<int, 2> state_idxs_solver = {
    0, 1
};

// clang-format on

} // namespace amici::model_model_dirac_py
