#include <amici/defines.h>
#include <array>
#include <string_view>

namespace amici::model_model_calvetti_py {

// clang-format off

extern const std::array<std::string_view const, 0> free_parameter_names = {
    
};

extern const std::array<std::string_view const, 6> fixed_parameter_names = {
    "V1ss", // k[0]
"R1ss", // k[1]
"V2ss", // k[2]
"R2ss", // k[3]
"V3ss", // k[4]
"R3ss", // k[5]
};

extern const std::array<std::string_view const, 6> state_names = {
    "V1", // x_rdata[0]
"V2", // x_rdata[1]
"V3", // x_rdata[2]
"f1", // x_rdata[3]
"f2", // x_rdata[4]
"f3", // x_rdata[5]
};

extern const std::array<std::string_view const, 6> state_names_solver = {
    "V1", // x_solver[0]
"V2", // x_solver[1]
"V3", // x_solver[2]
"f1", // x_solver[3]
"f2", // x_solver[4]
"f3", // x_solver[5]
};

extern const std::array<std::string_view const, 6> observable_names = {
    "y0", // y[0]
"y1", // y[1]
"y2", // y[2]
"y3", // y[3]
"y4", // y[4]
"y5", // y[5]
};

std::array<const ObservableScaling, 6> observable_scalings = {
    ObservableScaling::lin, // y[0]
ObservableScaling::lin, // y[1]
ObservableScaling::lin, // y[2]
ObservableScaling::lin, // y[3]
ObservableScaling::lin, // y[4]
ObservableScaling::lin, // y[5]
};

extern const std::array<std::string_view const, 16> expression_names = {
    "C1ss", // w[0]
"C2ss", // w[1]
"C3ss", // w[2]
"L1", // w[3]
"L2", // w[4]
"L3", // w[5]
"p2", // w[6]
"p3", // w[7]
"s", // w[8]
"R1", // w[9]
"R2", // w[10]
"R3", // w[11]
"f0", // w[12]
"rate_of_V1", // w[13]
"rate_of_V2", // w[14]
"rate_of_V3", // w[15]
};

extern const std::array<std::string_view const, 0> free_parameter_ids = {
    
};

extern const std::array<std::string_view const, 6> fixed_parameter_ids = {
    "V1ss", // k[0]
"R1ss", // k[1]
"V2ss", // k[2]
"R2ss", // k[3]
"V3ss", // k[4]
"R3ss", // k[5]
};

extern const std::array<std::string_view const, 6> state_ids = {
    "V1", // x_rdata[0]
"V2", // x_rdata[1]
"V3", // x_rdata[2]
"f1", // x_rdata[3]
"f2", // x_rdata[4]
"f3", // x_rdata[5]
};

extern const std::array<std::string_view const, 6> state_ids_solver = {
    "V1", // x_solver[0]
"V2", // x_solver[1]
"V3", // x_solver[2]
"f1", // x_solver[3]
"f2", // x_solver[4]
"f3", // x_solver[5]
};

extern const std::array<std::string_view const, 6> observable_ids = {
    "obs_V1", // y[0]
"obs_V2", // y[1]
"obs_V3", // y[2]
"obs_f0", // y[3]
"obs_f1", // y[4]
"obs_f2", // y[5]
};

extern const std::array<std::string_view const, 16> expression_ids = {
    "C1ss", // w[0]
"C2ss", // w[1]
"C3ss", // w[2]
"L1", // w[3]
"L2", // w[4]
"L3", // w[5]
"p2", // w[6]
"p3", // w[7]
"s", // w[8]
"R1", // w[9]
"R2", // w[10]
"R3", // w[11]
"f0", // w[12]
"rate_of_V1", // w[13]
"rate_of_V2", // w[14]
"rate_of_V3", // w[15]
};

std::array<int, 6> state_idxs_solver = {
    0, 1, 2, 3, 4, 5
};

// clang-format on

} // namespace amici::model_model_calvetti_py
