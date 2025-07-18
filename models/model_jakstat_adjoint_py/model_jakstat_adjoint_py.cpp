#include <amici/defines.h>
#include <array>

namespace amici {

namespace model_model_jakstat_adjoint_py {

// clang-format off

std::array<const char*, 17> parameterNames = {
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

std::array<const char*, 2> fixedParameterNames = {
    "Omega_cyt", // k[0]
"Omega_nuc", // k[1]
};

std::array<const char*, 9> stateNames = {
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

std::array<const char*, 3> observableNames = {
    "y0", // y[0]
"y1", // y[1]
"y2", // y[2]
};

std::array<const ObservableScaling, 3> observableScalings = {
    ObservableScaling::lin, // y[0]
ObservableScaling::lin, // y[1]
ObservableScaling::lin, // y[2]
};

std::array<const char*, 2> expressionNames = {
    "u", // w[0]
"flux_r0", // w[1]
};

std::array<const char*, 17> parameterIds = {
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

std::array<const char*, 2> fixedParameterIds = {
    "Omega_cyt", // k[0]
"Omega_nuc", // k[1]
};

std::array<const char*, 9> stateIds = {
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

std::array<const char*, 3> observableIds = {
    "obs_pSTAT", // y[0]
"obs_tSTAT", // y[1]
"obs_spline", // y[2]
};

std::array<const char*, 2> expressionIds = {
    "u", // w[0]
"flux_r0", // w[1]
};

std::array<int, 9> stateIdxsSolver = {
    0, 1, 2, 3, 4, 5, 6, 7, 8
};

// clang-format on

} // namespace model_model_jakstat_adjoint_py

} // namespace amici
