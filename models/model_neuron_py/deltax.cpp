#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"
#include "xdot.h"
#include "x_old.h"

namespace amici {
namespace model_model_neuron_py {

void deltax_model_neuron_py(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old){
    switch(ie) {
        case 0:
            deltax[0] = -c - v;
            deltax[1] = d - u + x_old1;
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
