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
#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"
#include "dxdotdp_explicit.h"

namespace amici {
namespace model_model_jakstat_adjoint_py {

void dxdotdp_explicit_model_jakstat_adjoint_py(realtype *dxdotdp_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddSTATdt_dp1 = -STAT*u;  // dxdotdp_explicit[0]
    ddpSTATdt_dp1 = STAT*u;  // dxdotdp_explicit[1]
    ddpSTATdt_dp2 = -2*std::pow(pSTAT, 2);  // dxdotdp_explicit[2]
    ddpSTAT_pSTATdt_dp2 = std::pow(pSTAT, 2);  // dxdotdp_explicit[3]
    ddpSTAT_pSTATdt_dp3 = -pSTAT_pSTAT;  // dxdotdp_explicit[4]
    ddnpSTAT_npSTATdt_dp3 = Omega_cyt*pSTAT_pSTAT/Omega_nuc;  // dxdotdp_explicit[5]
    ddSTATdt_dp4 = Omega_nuc*nSTAT5/Omega_cyt;  // dxdotdp_explicit[6]
    ddnpSTAT_npSTATdt_dp4 = -npSTAT_npSTAT;  // dxdotdp_explicit[7]
    ddnSTAT1dt_dp4 = -nSTAT1 + 2*npSTAT_npSTAT;  // dxdotdp_explicit[8]
    ddnSTAT2dt_dp4 = nSTAT1 - nSTAT2;  // dxdotdp_explicit[9]
    ddnSTAT3dt_dp4 = nSTAT2 - nSTAT3;  // dxdotdp_explicit[10]
    ddnSTAT4dt_dp4 = nSTAT3 - nSTAT4;  // dxdotdp_explicit[11]
    ddnSTAT5dt_dp4 = nSTAT4 - nSTAT5;  // dxdotdp_explicit[12]
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
