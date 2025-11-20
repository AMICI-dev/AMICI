#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_robertson_py {

static constexpr std::array<sunindextype, 4> dxdotdx_explicit_colptrs_model_robertson_py_ = {
    0, 3, 6, 9
};

void dxdotdx_explicit_colptrs_model_robertson_py(SUNMatrixWrapper &dxdotdx_explicit){
    dxdotdx_explicit.set_indexptrs(gsl::make_span(dxdotdx_explicit_colptrs_model_robertson_py_));
}
} // namespace model_model_robertson_py
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_robertson_py {

static constexpr std::array<sunindextype, 9> dxdotdx_explicit_rowvals_model_robertson_py_ = {
    0, 1, 2, 0, 1, 2, 0, 1, 2
};

void dxdotdx_explicit_rowvals_model_robertson_py(SUNMatrixWrapper &dxdotdx_explicit){
    dxdotdx_explicit.set_indexvals(gsl::make_span(dxdotdx_explicit_rowvals_model_robertson_py_));
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
#include "dxdotdx_explicit.h"

namespace amici {
namespace model_model_robertson_py {

void dxdotdx_explicit_model_robertson_py(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *w){
    dde_0_dx1 = -p1;  // dxdotdx_explicit[0]
    dde_1_dx1 = p1;  // dxdotdx_explicit[1]
    dae_0_dx1 = 1;  // dxdotdx_explicit[2]
    dde_0_dx2 = p2*x3;  // dxdotdx_explicit[3]
    dde_1_dx2 = -p2*x3 - 2*p3*x2;  // dxdotdx_explicit[4]
    dae_0_dx2 = 1;  // dxdotdx_explicit[5]
    dde_0_dx3 = p2*x2;  // dxdotdx_explicit[6]
    dde_1_dx3 = -p2*x2;  // dxdotdx_explicit[7]
    dae_0_dx3 = 1;  // dxdotdx_explicit[8]
}

} // namespace model_model_robertson_py
} // namespace amici
