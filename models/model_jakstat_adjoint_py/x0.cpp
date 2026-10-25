#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_jakstat_adjoint_py {

void x0_model_jakstat_adjoint_py(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    const realtype init_STAT_ = p[4];

    x0[0] = init_STAT_;
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
