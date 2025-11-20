#include <amici/defines.h>
#include <array>

namespace amici::model_model_jakstat_adjoint_py {

// clang-format off

std::array<const char*, 17> free_parameter_names = {
    "p1", // p[0]
"p2", // p[1]
"p3", // p[2]
"p4", // p[3]
"init_STAT", // p[4]
"sp1", // p[5]
"sp2", // p[6]
"sp3", // p[7]
"sp4", // p[8]
"sp5", // p[9]
"offset_tSTAT", // p[10]
"offset_pSTAT", // p[11]
"scale_tSTAT", // p[12]
"scale_pSTAT", // p[13]
"sigma_pSTAT", // p[14]
"sigma_tSTAT", // p[15]
"sigma_pEpoR", // p[16]
};

std::array<const char*, 2> fixed_parameter_names = {
    "Omega_cyt", // k[0]
"Omega_nuc", // k[1]
};

std::array<const char*, 9> state_names = {
    "STAT", // x_rdata[0]
"pSTAT", // x_rdata[1]
"pSTAT_pSTAT", // x_rdata[2]
"npSTAT_npSTAT", // x_rdata[3]
"nSTAT1", // x_rdata[4]
"nSTAT2", // x_rdata[5]
"nSTAT3", // x_rdata[6]
"nSTAT4", // x_rdata[7]
"nSTAT5", // x_rdata[8]
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

std::array<const char*, 1> expression_names = {
    "u", // w[0]
};

std::array<const char*, 17> free_parameter_ids = {
    "p1", // p[0]
"p2", // p[1]
"p3", // p[2]
"p4", // p[3]
"init_STAT", // p[4]
"sp1", // p[5]
"sp2", // p[6]
"sp3", // p[7]
"sp4", // p[8]
"sp5", // p[9]
"offset_tSTAT", // p[10]
"offset_pSTAT", // p[11]
"scale_tSTAT", // p[12]
"scale_pSTAT", // p[13]
"sigma_pSTAT", // p[14]
"sigma_tSTAT", // p[15]
"sigma_pEpoR", // p[16]
};

std::array<const char*, 2> fixed_parameter_ids = {
    "Omega_cyt", // k[0]
"Omega_nuc", // k[1]
};

std::array<const char*, 9> state_ids = {
    "STAT", // x_rdata[0]
"pSTAT", // x_rdata[1]
"pSTAT_pSTAT", // x_rdata[2]
"npSTAT_npSTAT", // x_rdata[3]
"nSTAT1", // x_rdata[4]
"nSTAT2", // x_rdata[5]
"nSTAT3", // x_rdata[6]
"nSTAT4", // x_rdata[7]
"nSTAT5", // x_rdata[8]
};

std::array<const char*, 3> observable_ids = {
    "obs_pSTAT", // y[0]
"obs_tSTAT", // y[1]
"obs_spline", // y[2]
};

std::array<const char*, 1> expression_ids = {
    "u", // w[0]
};

std::array<int, 9> state_idxs_solver = {
    0, 1, 2, 3, 4, 5, 6, 7, 8
};

// clang-format on

} // namespace amici::model_model_jakstat_adjoint_py
