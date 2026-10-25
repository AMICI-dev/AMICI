#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_dirac_py {

static constexpr std::array<sunindextype, 3> dxdotdx_explicit_colptrs_model_dirac_py_ = {
    0, 2, 3
};

void dxdotdx_explicit_colptrs_model_dirac_py(SUNMatrixWrapper &dxdotdx_explicit){
    dxdotdx_explicit.set_indexptrs(gsl::make_span(dxdotdx_explicit_colptrs_model_dirac_py_));
}
} // namespace model_model_dirac_py
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_dirac_py {

static constexpr std::array<sunindextype, 3> dxdotdx_explicit_rowvals_model_dirac_py_ = {
    0, 1, 1
};

void dxdotdx_explicit_rowvals_model_dirac_py(SUNMatrixWrapper &dxdotdx_explicit){
    dxdotdx_explicit.set_indexvals(gsl::make_span(dxdotdx_explicit_rowvals_model_dirac_py_));
}
} // namespace model_model_dirac_py
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include <sundials/sundials_types.h>
#include <gsl/gsl-lite.hpp>

namespace amici {
namespace model_model_dirac_py {

void dxdotdx_explicit_model_dirac_py(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype p1_ = p[0];
    const realtype p3_ = p[2];
    const realtype p4_ = p[3];

    realtype &ddx1dt_dx1_ = dxdotdx_explicit[0];
    realtype &ddx2dt_dx1_ = dxdotdx_explicit[1];
    realtype &ddx2dt_dx2_ = dxdotdx_explicit[2];
    ddx1dt_dx1_ = -p1_;  // dxdotdx_explicit[0]
    ddx2dt_dx1_ = p3_;  // dxdotdx_explicit[1]
    ddx2dt_dx2_ = -p4_;  // dxdotdx_explicit[2]
}

} // namespace model_model_dirac_py
} // namespace amici
