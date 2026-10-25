#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_nested_events_py {

static constexpr std::array<sunindextype, 6> dxdotdp_explicit_colptrs_model_nested_events_py_ = {
    0, 0, 0, 0, 1, 2
};

void dxdotdp_explicit_colptrs_model_nested_events_py(SUNMatrixWrapper &dxdotdp_explicit){
    dxdotdp_explicit.set_indexptrs(gsl::make_span(dxdotdp_explicit_colptrs_model_nested_events_py_));
}
} // namespace model_model_nested_events_py
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_nested_events_py {

static constexpr std::array<sunindextype, 2> dxdotdp_explicit_rowvals_model_nested_events_py_ = {
    0, 0
};

void dxdotdp_explicit_rowvals_model_nested_events_py(SUNMatrixWrapper &dxdotdp_explicit){
    dxdotdp_explicit.set_indexvals(gsl::make_span(dxdotdp_explicit_rowvals_model_nested_events_py_));
}
} // namespace model_model_nested_events_py
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include <sundials/sundials_types.h>
#include <gsl/gsl-lite.hpp>

namespace amici {
namespace model_model_nested_events_py {

void dxdotdp_explicit_model_nested_events_py(realtype *dxdotdp_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype Virus_ = x[0];
    const realtype Heaviside_1_ = h[0];

    realtype &ddVirusdt_drho_V_ = dxdotdp_explicit[0];
    realtype &ddVirusdt_ddelta_V_ = dxdotdp_explicit[1];
    ddVirusdt_drho_V_ = Heaviside_1_*Virus_;  // dxdotdp_explicit[0]
    ddVirusdt_ddelta_V_ = -Virus_;  // dxdotdp_explicit[1]
}

} // namespace model_model_nested_events_py
} // namespace amici
