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
#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"
#include "dxdotdp_explicit.h"

namespace amici {
namespace model_model_events_py {

void dxdotdp_explicit_model_events_py(realtype *dxdotdp_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddx1dt_dp1 = -x1*(1 - Heaviside_2);  // dxdotdp_explicit[0]
    ddx2dt_dp2 = x1*std::exp(-1.0/10.0*t);  // dxdotdp_explicit[1]
    ddx2dt_dp3 = -x2;  // dxdotdp_explicit[2]
}

} // namespace model_model_events_py
} // namespace amici
