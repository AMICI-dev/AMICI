#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_jakstat_adjoint_py {

static constexpr std::array<sunindextype, 18> dxdotdp_explicit_colptrs_model_jakstat_adjoint_py_ = {
    0, 2, 4, 6, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13
};

void dxdotdp_explicit_colptrs_model_jakstat_adjoint_py(SUNMatrixWrapper &dxdotdp_explicit){
    dxdotdp_explicit.set_indexptrs(gsl::make_span(dxdotdp_explicit_colptrs_model_jakstat_adjoint_py_));
}
} // namespace model_model_jakstat_adjoint_py
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_jakstat_adjoint_py {

static constexpr std::array<sunindextype, 13> dxdotdp_explicit_rowvals_model_jakstat_adjoint_py_ = {
    0, 1, 1, 2, 2, 3, 0, 3, 4, 5, 6, 7, 8
};

void dxdotdp_explicit_rowvals_model_jakstat_adjoint_py(SUNMatrixWrapper &dxdotdp_explicit){
    dxdotdp_explicit.set_indexvals(gsl::make_span(dxdotdp_explicit_rowvals_model_jakstat_adjoint_py_));
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

void dxdotdp_explicit_model_jakstat_adjoint_py(realtype *dxdotdp_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype STAT_ = x[0];
    const realtype pSTAT_ = x[1];
    const realtype pSTAT_pSTAT_ = x[2];
    const realtype npSTAT_npSTAT_ = x[3];
    const realtype nSTAT1_ = x[4];
    const realtype nSTAT2_ = x[5];
    const realtype nSTAT3_ = x[6];
    const realtype nSTAT4_ = x[7];
    const realtype nSTAT5_ = x[8];
    const realtype Omega_cyt_ = k[0];
    const realtype Omega_nuc_ = k[1];
    const realtype u_ = w[0];

    realtype &ddSTATdt_dp1_ = dxdotdp_explicit[0];
    realtype &ddpSTATdt_dp1_ = dxdotdp_explicit[1];
    realtype &ddpSTATdt_dp2_ = dxdotdp_explicit[2];
    realtype &ddpSTAT_pSTATdt_dp2_ = dxdotdp_explicit[3];
    realtype &ddpSTAT_pSTATdt_dp3_ = dxdotdp_explicit[4];
    realtype &ddnpSTAT_npSTATdt_dp3_ = dxdotdp_explicit[5];
    realtype &ddSTATdt_dp4_ = dxdotdp_explicit[6];
    realtype &ddnpSTAT_npSTATdt_dp4_ = dxdotdp_explicit[7];
    realtype &ddnSTAT1dt_dp4_ = dxdotdp_explicit[8];
    realtype &ddnSTAT2dt_dp4_ = dxdotdp_explicit[9];
    realtype &ddnSTAT3dt_dp4_ = dxdotdp_explicit[10];
    realtype &ddnSTAT4dt_dp4_ = dxdotdp_explicit[11];
    realtype &ddnSTAT5dt_dp4_ = dxdotdp_explicit[12];
    ddSTATdt_dp1_ = -STAT_*u_;  // dxdotdp_explicit[0]
    ddpSTATdt_dp1_ = STAT_*u_;  // dxdotdp_explicit[1]
    ddpSTATdt_dp2_ = -2*std::pow(pSTAT_, 2);  // dxdotdp_explicit[2]
    ddpSTAT_pSTATdt_dp2_ = std::pow(pSTAT_, 2);  // dxdotdp_explicit[3]
    ddpSTAT_pSTATdt_dp3_ = -pSTAT_pSTAT_;  // dxdotdp_explicit[4]
    ddnpSTAT_npSTATdt_dp3_ = Omega_cyt_*pSTAT_pSTAT_/Omega_nuc_;  // dxdotdp_explicit[5]
    ddSTATdt_dp4_ = Omega_nuc_*nSTAT5_/Omega_cyt_;  // dxdotdp_explicit[6]
    ddnpSTAT_npSTATdt_dp4_ = -npSTAT_npSTAT_;  // dxdotdp_explicit[7]
    ddnSTAT1dt_dp4_ = -nSTAT1_ + 2*npSTAT_npSTAT_;  // dxdotdp_explicit[8]
    ddnSTAT2dt_dp4_ = nSTAT1_ - nSTAT2_;  // dxdotdp_explicit[9]
    ddnSTAT3dt_dp4_ = nSTAT2_ - nSTAT3_;  // dxdotdp_explicit[10]
    ddnSTAT4dt_dp4_ = nSTAT3_ - nSTAT4_;  // dxdotdp_explicit[11]
    ddnSTAT5dt_dp4_ = nSTAT4_ - nSTAT5_;  // dxdotdp_explicit[12]
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
