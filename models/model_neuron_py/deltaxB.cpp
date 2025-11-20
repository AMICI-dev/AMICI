#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"
#include "xdot.h"
#include "xdot_old.h"
#include "x_old.h"
#include "xB.h"

namespace amici {
namespace model_model_neuron_py {

void deltaxB_model_neuron_py(realtype *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old, const realtype *xB, const realtype *tcl){
    switch(ie) {
        case 0:
            deltaxB[0] = xB0*(xdot_old0/(I0 - u + (1.0/25.0)*std::pow(v, 2) + 5*v + 140) + (dvdt - xdot_old0)/(I0 - u + (1.0/25.0)*std::pow(v, 2) + 5*v + 140) - 1) + xB1*(dudt - xdot_old1)/(I0 - u + (1.0/25.0)*std::pow(v, 2) + 5*v + 140);
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
