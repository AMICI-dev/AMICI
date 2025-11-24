#include <amici/defines.h>
#include <array>

namespace amici::model_model_neuron_py {

// clang-format off

std::array<const char*, 4> free_parameter_names = {
    "a", // p[0]
"b", // p[1]
"c", // p[2]
"d", // p[3]
};

std::array<const char*, 2> fixed_parameter_names = {
    "v0", // k[0]
"I0", // k[1]
};

std::array<const char*, 2> state_names = {
    "v", // x_rdata[0]
"u", // x_rdata[1]
};

std::array<const char*, 1> observable_names = {
    "v", // y[0]
};

std::array<const ObservableScaling, 1> observable_scalings = {
    ObservableScaling::lin, // y[0]
};

std::array<const char*, 0> expression_names = {
    
};

std::array<const char*, 4> free_parameter_ids = {
    "a", // p[0]
"b", // p[1]
"c", // p[2]
"d", // p[3]
};

std::array<const char*, 2> fixed_parameter_ids = {
    "v0", // k[0]
"I0", // k[1]
};

std::array<const char*, 2> state_ids = {
    "v", // x_rdata[0]
"u", // x_rdata[1]
};

std::array<const char*, 1> observable_ids = {
    "y1", // y[0]
};

std::array<const char*, 0> expression_ids = {
    
};

std::array<int, 2> state_idxs_solver = {
    0, 1
};

// clang-format on

} // namespace amici::model_model_neuron_py
