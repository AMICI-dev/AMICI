#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_events_py {

static constexpr std::array<sunindextype, 5> dxdotdp_explicit_colptrs_model_events_py_ = {
    0, 1, 2, 3, 3
};

void dxdotdp_explicit_colptrs_model_events_py(SUNMatrixWrapper &dxdotdp_explicit){
    dxdotdp_explicit.set_indexptrs(gsl::make_span(dxdotdp_explicit_colptrs_model_events_py_));
}
} // namespace model_model_events_py
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_events_py {

static constexpr std::array<sunindextype, 3> dxdotdp_explicit_rowvals_model_events_py_ = {
    0, 1, 1
};

void dxdotdp_explicit_rowvals_model_events_py(SUNMatrixWrapper &dxdotdp_explicit){
    dxdotdp_explicit.set_indexvals(gsl::make_span(dxdotdp_explicit_rowvals_model_events_py_));
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

void dxdotdp_explicit_model_events_py(realtype *dxdotdp_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype x1_ = x[0];
    const realtype x2_ = x[1];
    const realtype Heaviside_2_ = h[2];

    realtype &ddx1dt_dp1_ = dxdotdp_explicit[0];
    realtype &ddx2dt_dp2_ = dxdotdp_explicit[1];
    realtype &ddx2dt_dp3_ = dxdotdp_explicit[2];
    ddx1dt_dp1_ = -x1_*(1 - Heaviside_2_);  // dxdotdp_explicit[0]
    ddx2dt_dp2_ = x1_*std::exp(-1.0/10.0*t);  // dxdotdp_explicit[1]
    ddx2dt_dp3_ = -x2_;  // dxdotdp_explicit[2]
}

} // namespace model_model_events_py
} // namespace amici
