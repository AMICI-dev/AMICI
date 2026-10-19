#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

void dtaudx_model_neuron_py(realtype *dtaudx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const realtype *tcl, const int ie){
    const realtype v_ = x[0];
    const realtype u_ = x[1];
    const realtype I0_ = k[1];

    switch(ie) {
        case 0:
            dtaudx[0] = 1.0/(I0_ - u_ + (1.0/25.0)*std::pow(v_, 2) + 5*v_ + 140);
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
