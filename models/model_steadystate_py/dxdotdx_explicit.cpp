#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_steadystate_py {

static constexpr std::array<sunindextype, 4> dxdotdx_explicit_colptrs_model_steadystate_py_ = {
    0, 3, 6, 9
};

void dxdotdx_explicit_colptrs_model_steadystate_py(SUNMatrixWrapper &dxdotdx_explicit){
    dxdotdx_explicit.set_indexptrs(gsl::make_span(dxdotdx_explicit_colptrs_model_steadystate_py_));
}
} // namespace model_model_steadystate_py
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_steadystate_py {

static constexpr std::array<sunindextype, 9> dxdotdx_explicit_rowvals_model_steadystate_py_ = {
    0, 1, 2, 0, 1, 2, 0, 1, 2
};

void dxdotdx_explicit_rowvals_model_steadystate_py(SUNMatrixWrapper &dxdotdx_explicit){
    dxdotdx_explicit.set_indexvals(gsl::make_span(dxdotdx_explicit_rowvals_model_steadystate_py_));
}
} // namespace model_model_steadystate_py
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include <sundials/sundials_types.h>
#include <gsl/gsl-lite.hpp>

namespace amici {
namespace model_model_steadystate_py {

void dxdotdx_explicit_model_steadystate_py(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype x1_ = x[0];
    const realtype x2_ = x[1];
    const realtype p1_ = p[0];
    const realtype p2_ = p[1];
    const realtype p3_ = p[2];
    const realtype p4_ = p[3];
    const realtype k4_ = k[3];

    realtype &ddx1dt_dx1_ = dxdotdx_explicit[0];
    realtype &ddx2dt_dx1_ = dxdotdx_explicit[1];
    realtype &ddx3dt_dx1_ = dxdotdx_explicit[2];
    realtype &ddx1dt_dx2_ = dxdotdx_explicit[3];
    realtype &ddx2dt_dx2_ = dxdotdx_explicit[4];
    realtype &ddx3dt_dx2_ = dxdotdx_explicit[5];
    realtype &ddx1dt_dx3_ = dxdotdx_explicit[6];
    realtype &ddx2dt_dx3_ = dxdotdx_explicit[7];
    realtype &ddx3dt_dx3_ = dxdotdx_explicit[8];
    ddx1dt_dx1_ = -4*p1_*x1_ - p2_*x2_;  // dxdotdx_explicit[0]
    ddx2dt_dx1_ = 2*p1_*x1_ - p2_*x2_;  // dxdotdx_explicit[1]
    ddx3dt_dx1_ = p2_*x2_;  // dxdotdx_explicit[2]
    ddx1dt_dx2_ = -p2_*x1_ + 2*p3_;  // dxdotdx_explicit[3]
    ddx2dt_dx2_ = -p2_*x1_ - p3_;  // dxdotdx_explicit[4]
    ddx3dt_dx2_ = p2_*x1_;  // dxdotdx_explicit[5]
    ddx1dt_dx3_ = p4_;  // dxdotdx_explicit[6]
    ddx2dt_dx3_ = p4_;  // dxdotdx_explicit[7]
    ddx3dt_dx3_ = -k4_ - p4_;  // dxdotdx_explicit[8]
}

} // namespace model_model_steadystate_py
} // namespace amici
