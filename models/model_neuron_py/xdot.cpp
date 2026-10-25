#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

void xdot_model_neuron_py(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype v_ = x[0];
    const realtype u_ = x[1];
    const realtype a_ = p[0];
    const realtype b_ = p[1];
    const realtype I0_ = k[1];

    realtype &dvdt_ = xdot[0];
    realtype &dudt_ = xdot[1];
    dvdt_ = I0_ - u_ + (1.0/25.0)*std::pow(v_, 2) + 5*v_ + 140;  // xdot[0]
    dudt_ = a_*(b_*v_ - u_);  // xdot[1]
}

} // namespace model_model_neuron_py
} // namespace amici
