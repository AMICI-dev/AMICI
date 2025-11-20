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
#include "x.h"
#include "p.h"
#include "k.h"
#include "dxdotdx_explicit.h"

namespace amici {
namespace model_model_steadystate_py {

void dxdotdx_explicit_model_steadystate_py(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddx1dt_dx1 = -4*p1*x1 - p2*x2;  // dxdotdx_explicit[0]
    ddx2dt_dx1 = 2*p1*x1 - p2*x2;  // dxdotdx_explicit[1]
    ddx3dt_dx1 = p2*x2;  // dxdotdx_explicit[2]
    ddx1dt_dx2 = -p2*x1 + 2*p3;  // dxdotdx_explicit[3]
    ddx2dt_dx2 = -p2*x1 - p3;  // dxdotdx_explicit[4]
    ddx3dt_dx2 = p2*x1;  // dxdotdx_explicit[5]
    ddx1dt_dx3 = p4;  // dxdotdx_explicit[6]
    ddx2dt_dx3 = p4;  // dxdotdx_explicit[7]
    ddx3dt_dx3 = -k4 - p4;  // dxdotdx_explicit[8]
}

} // namespace model_model_steadystate_py
} // namespace amici
