#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_robertson_py {

static constexpr std::array<sunindextype, 4> dxdotdp_explicit_colptrs_model_robertson_py_ = {
    0, 2, 4, 5
};

void dxdotdp_explicit_colptrs_model_robertson_py(SUNMatrixWrapper &dxdotdp_explicit){
    dxdotdp_explicit.set_indexptrs(gsl::make_span(dxdotdp_explicit_colptrs_model_robertson_py_));
}
} // namespace model_model_robertson_py
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_robertson_py {

static constexpr std::array<sunindextype, 5> dxdotdp_explicit_rowvals_model_robertson_py_ = {
    0, 1, 0, 1, 1
};

void dxdotdp_explicit_rowvals_model_robertson_py(SUNMatrixWrapper &dxdotdp_explicit){
    dxdotdp_explicit.set_indexvals(gsl::make_span(dxdotdp_explicit_rowvals_model_robertson_py_));
}
} // namespace model_model_robertson_py
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include <sundials/sundials_types.h>
#include <gsl/gsl-lite.hpp>

namespace amici {
namespace model_model_robertson_py {

void dxdotdp_explicit_model_robertson_py(realtype *dxdotdp_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *w){
    const realtype x1_ = x[0];
    const realtype x2_ = x[1];
    const realtype x3_ = x[2];

    realtype &dde_0_dp1_ = dxdotdp_explicit[0];
    realtype &dde_1_dp1_ = dxdotdp_explicit[1];
    realtype &dde_0_dp2_ = dxdotdp_explicit[2];
    realtype &dde_1_dp2_ = dxdotdp_explicit[3];
    realtype &dde_1_dp3_ = dxdotdp_explicit[4];
    dde_0_dp1_ = -x1_;  // dxdotdp_explicit[0]
    dde_1_dp1_ = x1_;  // dxdotdp_explicit[1]
    dde_0_dp2_ = x2_*x3_;  // dxdotdp_explicit[2]
    dde_1_dp2_ = -x2_*x3_;  // dxdotdp_explicit[3]
    dde_1_dp3_ = -std::pow(x2_, 2);  // dxdotdp_explicit[4]
}

} // namespace model_model_robertson_py
} // namespace amici
