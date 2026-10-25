#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_jakstat_adjoint_py {

static constexpr std::array<sunindextype, 10> dxdotdx_explicit_colptrs_model_jakstat_adjoint_py_ = {
    0, 2, 4, 6, 8, 10, 12, 14, 16, 18
};

void dxdotdx_explicit_colptrs_model_jakstat_adjoint_py(SUNMatrixWrapper &dxdotdx_explicit){
    dxdotdx_explicit.set_indexptrs(gsl::make_span(dxdotdx_explicit_colptrs_model_jakstat_adjoint_py_));
}
} // namespace model_model_jakstat_adjoint_py
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_jakstat_adjoint_py {

static constexpr std::array<sunindextype, 18> dxdotdx_explicit_rowvals_model_jakstat_adjoint_py_ = {
    0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 0, 8
};

void dxdotdx_explicit_rowvals_model_jakstat_adjoint_py(SUNMatrixWrapper &dxdotdx_explicit){
    dxdotdx_explicit.set_indexvals(gsl::make_span(dxdotdx_explicit_rowvals_model_jakstat_adjoint_py_));
}
} // namespace model_model_jakstat_adjoint_py
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include <sundials/sundials_types.h>
#include <gsl/gsl-lite.hpp>

namespace amici {
namespace model_model_jakstat_adjoint_py {

void dxdotdx_explicit_model_jakstat_adjoint_py(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype pSTAT_ = x[1];
    const realtype p1_ = p[0];
    const realtype p2_ = p[1];
    const realtype p3_ = p[2];
    const realtype p4_ = p[3];
    const realtype Omega_cyt_ = k[0];
    const realtype Omega_nuc_ = k[1];
    const realtype u_ = w[0];

    realtype &ddSTATdt_dSTAT_ = dxdotdx_explicit[0];
    realtype &ddpSTATdt_dSTAT_ = dxdotdx_explicit[1];
    realtype &ddpSTATdt_dpSTAT_ = dxdotdx_explicit[2];
    realtype &ddpSTAT_pSTATdt_dpSTAT_ = dxdotdx_explicit[3];
    realtype &ddpSTAT_pSTATdt_dpSTAT_pSTAT_ = dxdotdx_explicit[4];
    realtype &ddnpSTAT_npSTATdt_dpSTAT_pSTAT_ = dxdotdx_explicit[5];
    realtype &ddnpSTAT_npSTATdt_dnpSTAT_npSTAT_ = dxdotdx_explicit[6];
    realtype &ddnSTAT1dt_dnpSTAT_npSTAT_ = dxdotdx_explicit[7];
    realtype &ddnSTAT1dt_dnSTAT1_ = dxdotdx_explicit[8];
    realtype &ddnSTAT2dt_dnSTAT1_ = dxdotdx_explicit[9];
    realtype &ddnSTAT2dt_dnSTAT2_ = dxdotdx_explicit[10];
    realtype &ddnSTAT3dt_dnSTAT2_ = dxdotdx_explicit[11];
    realtype &ddnSTAT3dt_dnSTAT3_ = dxdotdx_explicit[12];
    realtype &ddnSTAT4dt_dnSTAT3_ = dxdotdx_explicit[13];
    realtype &ddnSTAT4dt_dnSTAT4_ = dxdotdx_explicit[14];
    realtype &ddnSTAT5dt_dnSTAT4_ = dxdotdx_explicit[15];
    realtype &ddSTATdt_dnSTAT5_ = dxdotdx_explicit[16];
    realtype &ddnSTAT5dt_dnSTAT5_ = dxdotdx_explicit[17];
    ddSTATdt_dSTAT_ = -p1_*u_;  // dxdotdx_explicit[0]
    ddpSTATdt_dSTAT_ = p1_*u_;  // dxdotdx_explicit[1]
    ddpSTATdt_dpSTAT_ = -4*p2_*pSTAT_;  // dxdotdx_explicit[2]
    ddpSTAT_pSTATdt_dpSTAT_ = 2*p2_*pSTAT_;  // dxdotdx_explicit[3]
    ddpSTAT_pSTATdt_dpSTAT_pSTAT_ = -p3_;  // dxdotdx_explicit[4]
    ddnpSTAT_npSTATdt_dpSTAT_pSTAT_ = Omega_cyt_*p3_/Omega_nuc_;  // dxdotdx_explicit[5]
    ddnpSTAT_npSTATdt_dnpSTAT_npSTAT_ = -p4_;  // dxdotdx_explicit[6]
    ddnSTAT1dt_dnpSTAT_npSTAT_ = 2*p4_;  // dxdotdx_explicit[7]
    ddnSTAT1dt_dnSTAT1_ = -p4_;  // dxdotdx_explicit[8]
    ddnSTAT2dt_dnSTAT1_ = p4_;  // dxdotdx_explicit[9]
    ddnSTAT2dt_dnSTAT2_ = -p4_;  // dxdotdx_explicit[10]
    ddnSTAT3dt_dnSTAT2_ = p4_;  // dxdotdx_explicit[11]
    ddnSTAT3dt_dnSTAT3_ = -p4_;  // dxdotdx_explicit[12]
    ddnSTAT4dt_dnSTAT3_ = p4_;  // dxdotdx_explicit[13]
    ddnSTAT4dt_dnSTAT4_ = -p4_;  // dxdotdx_explicit[14]
    ddnSTAT5dt_dnSTAT4_ = p4_;  // dxdotdx_explicit[15]
    ddSTATdt_dnSTAT5_ = Omega_nuc_*p4_/Omega_cyt_;  // dxdotdx_explicit[16]
    ddnSTAT5dt_dnSTAT5_ = -p4_;  // dxdotdx_explicit[17]
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
