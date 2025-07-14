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
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"
#include "dxdotdx_explicit.h"

namespace amici {
namespace model_model_jakstat_adjoint_py {

void dxdotdx_explicit_model_jakstat_adjoint_py(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddSTATdt_dSTAT = -p1*u;  // dxdotdx_explicit[0]
    ddpSTATdt_dSTAT = p1*u;  // dxdotdx_explicit[1]
    ddpSTATdt_dpSTAT = -4*p2*pSTAT;  // dxdotdx_explicit[2]
    ddpSTAT_pSTATdt_dpSTAT = 2*p2*pSTAT;  // dxdotdx_explicit[3]
    ddpSTAT_pSTATdt_dpSTAT_pSTAT = -p3;  // dxdotdx_explicit[4]
    ddnpSTAT_npSTATdt_dpSTAT_pSTAT = Omega_cyt*p3/Omega_nuc;  // dxdotdx_explicit[5]
    ddnpSTAT_npSTATdt_dnpSTAT_npSTAT = -p4;  // dxdotdx_explicit[6]
    ddnSTAT1dt_dnpSTAT_npSTAT = 2*p4;  // dxdotdx_explicit[7]
    ddnSTAT1dt_dnSTAT1 = -p4;  // dxdotdx_explicit[8]
    ddnSTAT2dt_dnSTAT1 = p4;  // dxdotdx_explicit[9]
    ddnSTAT2dt_dnSTAT2 = -p4;  // dxdotdx_explicit[10]
    ddnSTAT3dt_dnSTAT2 = p4;  // dxdotdx_explicit[11]
    ddnSTAT3dt_dnSTAT3 = -p4;  // dxdotdx_explicit[12]
    ddnSTAT4dt_dnSTAT3 = p4;  // dxdotdx_explicit[13]
    ddnSTAT4dt_dnSTAT4 = -p4;  // dxdotdx_explicit[14]
    ddnSTAT5dt_dnSTAT4 = p4;  // dxdotdx_explicit[15]
    ddSTATdt_dnSTAT5 = Omega_nuc*p4/Omega_cyt;  // dxdotdx_explicit[16]
    ddnSTAT5dt_dnSTAT5 = -p4;  // dxdotdx_explicit[17]
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
