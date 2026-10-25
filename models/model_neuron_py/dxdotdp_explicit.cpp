#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

static constexpr std::array<sunindextype, 5> dxdotdp_explicit_colptrs_model_neuron_py_ = {
    0, 1, 2, 2, 2
};

void dxdotdp_explicit_colptrs_model_neuron_py(SUNMatrixWrapper &dxdotdp_explicit){
    dxdotdp_explicit.set_indexptrs(gsl::make_span(dxdotdp_explicit_colptrs_model_neuron_py_));
}
} // namespace model_model_neuron_py
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

static constexpr std::array<sunindextype, 2> dxdotdp_explicit_rowvals_model_neuron_py_ = {
    1, 1
};

void dxdotdp_explicit_rowvals_model_neuron_py(SUNMatrixWrapper &dxdotdp_explicit){
    dxdotdp_explicit.set_indexvals(gsl::make_span(dxdotdp_explicit_rowvals_model_neuron_py_));
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

void dxdotdp_explicit_model_neuron_py(realtype *dxdotdp_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype v_ = x[0];
    const realtype u_ = x[1];
    const realtype a_ = p[0];
    const realtype b_ = p[1];

    realtype &ddudt_da_ = dxdotdp_explicit[0];
    realtype &ddudt_db_ = dxdotdp_explicit[1];
    ddudt_da_ = b_*v_ - u_;  // dxdotdp_explicit[0]
    ddudt_db_ = a_*v_;  // dxdotdp_explicit[1]
}

} // namespace model_model_neuron_py
} // namespace amici
