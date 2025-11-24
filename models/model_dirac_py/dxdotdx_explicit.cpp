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
#include "x.h"
#include "p.h"
#include "h.h"
#include "dxdotdx_explicit.h"

namespace amici {
namespace model_model_dirac_py {

void dxdotdx_explicit_model_dirac_py(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddx1dt_dx1 = -p1;  // dxdotdx_explicit[0]
    ddx2dt_dx1 = p3;  // dxdotdx_explicit[1]
    ddx2dt_dx2 = -p4;  // dxdotdx_explicit[2]
}

} // namespace model_model_dirac_py
} // namespace amici
