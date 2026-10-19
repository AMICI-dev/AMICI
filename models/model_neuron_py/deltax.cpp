#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

void deltax_model_neuron_py(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_eval, const realtype *x_eval){
    const realtype v_ = x[0];
    const realtype u_ = x[1];
    const realtype c_ = p[2];
    const realtype d_ = p[3];
    const realtype x_eval1_ = x_eval[1];

    switch(ie) {
        case 0:
            deltax[0] = -c_ - v_;
            deltax[1] = d_ - u_ + x_eval1_;
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
