#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_robertson_py {

void x0_model_robertson_py(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    const realtype k1_ = k[0];

    x0[0] = k1_;
}

} // namespace model_model_robertson_py
} // namespace amici
