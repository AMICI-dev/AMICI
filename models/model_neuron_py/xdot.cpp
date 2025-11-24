#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"
#include "xdot.h"

namespace amici {
namespace model_model_neuron_py {

void xdot_model_neuron_py(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dvdt = I0 - u + (1.0/25.0)*std::pow(v, 2) + 5*v + 140;  // xdot[0]
    dudt = a*(b*v - u);  // xdot[1]
}

} // namespace model_model_neuron_py
} // namespace amici
