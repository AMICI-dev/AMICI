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
#include "x.h"
#include "p.h"
#include "k.h"
#include "dx.h"
#include "dxdotdp_explicit.h"

namespace amici {
namespace model_model_robertson_py {

void dxdotdp_explicit_model_robertson_py(realtype *dxdotdp_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *w){
    dde_0_dp1 = -x1;  // dxdotdp_explicit[0]
    dde_1_dp1 = x1;  // dxdotdp_explicit[1]
    dde_0_dp2 = x2*x3;  // dxdotdp_explicit[2]
    dde_1_dp2 = -x2*x3;  // dxdotdp_explicit[3]
    dde_1_dp3 = -std::pow(x2, 2);  // dxdotdp_explicit[4]
}

} // namespace model_model_robertson_py
} // namespace amici
