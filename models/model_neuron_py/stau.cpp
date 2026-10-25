#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

void stau_model_neuron_py(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const realtype *tcl, const realtype *sx, const int ip, const int ie){
    const realtype v_ = x[0];
    const realtype u_ = x[1];
    const realtype I0_ = k[1];
    const realtype sx0_ = sx[0];

    switch(ie) {
        case 0:
            switch(ip) {
                case 0:
                case 1:
                case 2:
                case 3:
                    stau[0] = sx0_/(I0_ - u_ + (1.0/25.0)*std::pow(v_, 2) + 5*v_ + 140);
                    break;
            }
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
