#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_jakstat_adjoint_py {

void w_model_jakstat_adjoint_py(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){
    const realtype spl_0_ = spl[0];

    realtype &u_ = w[0];

    // dynamic expressions
    u_ = spl_0_;  // w[0]
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
