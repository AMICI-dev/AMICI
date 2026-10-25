#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

static constexpr std::array<sunindextype, 3> dxdotdx_explicit_colptrs_model_neuron_py_ = {
    0, 2, 4
};

void dxdotdx_explicit_colptrs_model_neuron_py(SUNMatrixWrapper &dxdotdx_explicit){
    dxdotdx_explicit.set_indexptrs(gsl::make_span(dxdotdx_explicit_colptrs_model_neuron_py_));
}
} // namespace model_model_neuron_py
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

static constexpr std::array<sunindextype, 4> dxdotdx_explicit_rowvals_model_neuron_py_ = {
    0, 1, 0, 1
};

void dxdotdx_explicit_rowvals_model_neuron_py(SUNMatrixWrapper &dxdotdx_explicit){
    dxdotdx_explicit.set_indexvals(gsl::make_span(dxdotdx_explicit_rowvals_model_neuron_py_));
}
} // namespace model_model_neuron_py
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include <sundials/sundials_types.h>
#include <gsl/gsl-lite.hpp>

namespace amici {
namespace model_model_neuron_py {

void dxdotdx_explicit_model_neuron_py(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype v_ = x[0];
    const realtype a_ = p[0];
    const realtype b_ = p[1];

    realtype &ddvdt_dv_ = dxdotdx_explicit[0];
    realtype &ddudt_dv_ = dxdotdx_explicit[1];
    realtype &ddvdt_du_ = dxdotdx_explicit[2];
    realtype &ddudt_du_ = dxdotdx_explicit[3];
    ddvdt_dv_ = (2.0/25.0)*v_ + 5;  // dxdotdx_explicit[0]
    ddudt_dv_ = a_*b_;  // dxdotdx_explicit[1]
    ddvdt_du_ = -1;  // dxdotdx_explicit[2]
    ddudt_du_ = -a_;  // dxdotdx_explicit[3]
}

} // namespace model_model_neuron_py
} // namespace amici
