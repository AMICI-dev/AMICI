#include <amici/defines.h>
#include <array>

namespace amici::model_model_events_py {

// clang-format off

std::array<const char*, 4> free_parameter_names = {
    "p1", // p[0]
"p2", // p[1]
"p3", // p[2]
"p4", // p[3]
};

std::array<const char*, 4> fixed_parameter_names = {
    "k1", // k[0]
"k2", // k[1]
"k3", // k[2]
"k4", // k[3]
};

std::array<const char*, 3> state_names = {
    "x1", // x_rdata[0]
"x2", // x_rdata[1]
"x3", // x_rdata[2]
};

std::array<const char*, 1> observable_names = {
    "y1", // y[0]
};

std::array<const ObservableScaling, 1> observable_scalings = {
    ObservableScaling::lin, // y[0]
};

std::array<const char*, 0> expression_names = {
    
};

std::array<const char*, 4> free_parameter_ids = {
    "p1", // p[0]
"p2", // p[1]
"p3", // p[2]
"p4", // p[3]
};

std::array<const char*, 4> fixed_parameter_ids = {
    "k1", // k[0]
"k2", // k[1]
"k3", // k[2]
"k4", // k[3]
};

std::array<const char*, 3> state_ids = {
    "x1", // x_rdata[0]
"x2", // x_rdata[1]
"x3", // x_rdata[2]
};

std::array<const char*, 1> observable_ids = {
    "y1", // y[0]
};

std::array<const char*, 0> expression_ids = {
    
};

std::array<int, 3> state_idxs_solver = {
    0, 1, 2
};

// clang-format on

} // namespace amici::model_model_events_py
