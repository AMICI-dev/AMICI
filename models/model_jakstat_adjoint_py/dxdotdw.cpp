#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_jakstat_adjoint_py {

static constexpr std::array<sunindextype, 2> dxdotdw_colptrs_model_jakstat_adjoint_py_ = {
    0, 2
};

void dxdotdw_colptrs_model_jakstat_adjoint_py(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_model_jakstat_adjoint_py_));
}
} // namespace model_model_jakstat_adjoint_py
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_jakstat_adjoint_py {

static constexpr std::array<sunindextype, 2> dxdotdw_rowvals_model_jakstat_adjoint_py_ = {
    0, 1
};

void dxdotdw_rowvals_model_jakstat_adjoint_py(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_model_jakstat_adjoint_py_));
}
} // namespace model_model_jakstat_adjoint_py
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include <sundials/sundials_types.h>
#include <gsl/gsl-lite.hpp>
#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"
#include "dxdotdw.h"

namespace amici {
namespace model_model_jakstat_adjoint_py {

void dxdotdw_model_jakstat_adjoint_py(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddSTATdt_du = -STAT*p1;  // dxdotdw[0]
    ddpSTATdt_du = STAT*p1;  // dxdotdw[1]
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
