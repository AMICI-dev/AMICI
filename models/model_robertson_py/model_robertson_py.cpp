#include <amici/defines.h>
#include <array>

namespace amici::model_model_robertson_py {

// clang-format off

std::array<const char*, 3> free_parameter_names = {
    "p1", // p[0]
"p2", // p[1]
"p3", // p[2]
};

std::array<const char*, 1> fixed_parameter_names = {
    "k1", // k[0]
};

std::array<const char*, 3> state_names = {
    "x1", // x_rdata[0]
"x2", // x_rdata[1]
"x3", // x_rdata[2]
};

std::array<const char*, 3> observable_names = {
    "y0", // y[0]
"y1", // y[1]
"y2", // y[2]
};

std::array<const ObservableScaling, 3> observable_scalings = {
    ObservableScaling::lin, // y[0]
ObservableScaling::lin, // y[1]
ObservableScaling::lin, // y[2]
};

std::array<const char*, 0> expression_names = {
    
};

std::array<const char*, 3> free_parameter_ids = {
    "p1", // p[0]
"p2", // p[1]
"p3", // p[2]
};

std::array<const char*, 1> fixed_parameter_ids = {
    "k1", // k[0]
};

std::array<const char*, 3> state_ids = {
    "x1", // x_rdata[0]
"x2", // x_rdata[1]
"x3", // x_rdata[2]
};

std::array<const char*, 3> observable_ids = {
    "obs_x1", // y[0]
"obs_x2", // y[1]
"obs_x3", // y[2]
};

std::array<const char*, 0> expression_ids = {
    
};

std::array<int, 3> state_idxs_solver = {
    0, 1, 2
};

// clang-format on

} // namespace amici::model_model_robertson_py
