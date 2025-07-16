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
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"
#include "w.h"
#include "dxdotdp_explicit.h"

namespace amici {
namespace model_model_neuron_py {

void dxdotdp_explicit_model_neuron_py(realtype *dxdotdp_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddudt_da = b*v - u;  // dxdotdp_explicit[0]
    ddudt_db = a*v;  // dxdotdp_explicit[1]
}

} // namespace model_model_neuron_py
} // namespace amici
