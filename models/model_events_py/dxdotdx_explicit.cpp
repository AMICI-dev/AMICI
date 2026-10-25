#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_events_py {

static constexpr std::array<sunindextype, 4> dxdotdx_explicit_colptrs_model_events_py_ = {
    0, 2, 3, 4
};

void dxdotdx_explicit_colptrs_model_events_py(SUNMatrixWrapper &dxdotdx_explicit){
    dxdotdx_explicit.set_indexptrs(gsl::make_span(dxdotdx_explicit_colptrs_model_events_py_));
}
} // namespace model_model_events_py
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_events_py {

static constexpr std::array<sunindextype, 4> dxdotdx_explicit_rowvals_model_events_py_ = {
    0, 1, 1, 2
};

void dxdotdx_explicit_rowvals_model_events_py(SUNMatrixWrapper &dxdotdx_explicit){
    dxdotdx_explicit.set_indexvals(gsl::make_span(dxdotdx_explicit_rowvals_model_events_py_));
}
} // namespace model_model_events_py
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include <sundials/sundials_types.h>
#include <gsl/gsl-lite.hpp>

namespace amici {
namespace model_model_events_py {

void dxdotdx_explicit_model_events_py(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype p1_ = p[0];
    const realtype p2_ = p[1];
    const realtype p3_ = p[2];
    const realtype Heaviside_2_ = h[2];

    realtype &ddx1dt_dx1_ = dxdotdx_explicit[0];
    realtype &ddx2dt_dx1_ = dxdotdx_explicit[1];
    realtype &ddx2dt_dx2_ = dxdotdx_explicit[2];
    realtype &ddx3dt_dx3_ = dxdotdx_explicit[3];
    ddx1dt_dx1_ = -p1_*(1 - Heaviside_2_);  // dxdotdx_explicit[0]
    ddx2dt_dx1_ = p2_*std::exp(-1.0/10.0*t);  // dxdotdx_explicit[1]
    ddx2dt_dx2_ = -p3_;  // dxdotdx_explicit[2]
    ddx3dt_dx3_ = -1;  // dxdotdx_explicit[3]
}

} // namespace model_model_events_py
} // namespace amici
