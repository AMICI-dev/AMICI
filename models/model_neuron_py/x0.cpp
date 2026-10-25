#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

void x0_model_neuron_py(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    const realtype b_ = p[1];
    const realtype v0_ = k[0];

    x0[0] = v0_;
    x0[1] = b_*v0_;
}

} // namespace model_model_neuron_py
} // namespace amici
