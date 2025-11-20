#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_steadystate_py {

static constexpr std::array<sunindextype, 6> dxdotdp_explicit_colptrs_model_steadystate_py_ = {
    0, 2, 5, 7, 10, 11
};

void dxdotdp_explicit_colptrs_model_steadystate_py(SUNMatrixWrapper &dxdotdp_explicit){
    dxdotdp_explicit.set_indexptrs(gsl::make_span(dxdotdp_explicit_colptrs_model_steadystate_py_));
}
} // namespace model_model_steadystate_py
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_steadystate_py {

static constexpr std::array<sunindextype, 11> dxdotdp_explicit_rowvals_model_steadystate_py_ = {
    0, 1, 0, 1, 2, 0, 1, 0, 1, 2, 0
};

void dxdotdp_explicit_rowvals_model_steadystate_py(SUNMatrixWrapper &dxdotdp_explicit){
    dxdotdp_explicit.set_indexvals(gsl::make_span(dxdotdp_explicit_rowvals_model_steadystate_py_));
}
} // namespace model_model_steadystate_py
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include <sundials/sundials_types.h>
#include <gsl/gsl-lite.hpp>
#include "x.h"
#include "p.h"
#include "k.h"
#include "dxdotdp_explicit.h"

namespace amici {
namespace model_model_steadystate_py {

void dxdotdp_explicit_model_steadystate_py(realtype *dxdotdp_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddx1dt_dp1 = -2*std::pow(x1, 2);  // dxdotdp_explicit[0]
    ddx2dt_dp1 = std::pow(x1, 2);  // dxdotdp_explicit[1]
    ddx1dt_dp2 = -x1*x2;  // dxdotdp_explicit[2]
    ddx2dt_dp2 = -x1*x2;  // dxdotdp_explicit[3]
    ddx3dt_dp2 = x1*x2;  // dxdotdp_explicit[4]
    ddx1dt_dp3 = 2*x2;  // dxdotdp_explicit[5]
    ddx2dt_dp3 = -x2;  // dxdotdp_explicit[6]
    ddx1dt_dp4 = x3;  // dxdotdp_explicit[7]
    ddx2dt_dp4 = x3;  // dxdotdp_explicit[8]
    ddx3dt_dp4 = -x3;  // dxdotdp_explicit[9]
    ddx1dt_dp5 = 1;  // dxdotdp_explicit[10]
}

} // namespace model_model_steadystate_py
} // namespace amici
