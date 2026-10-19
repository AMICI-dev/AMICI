#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

void deltaqB_model_neuron_py(realtype *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const int ip, const int ie, const realtype *xdot, const realtype *xdot_eval, const realtype *x_eval, const realtype *xB){
    const realtype xB0_ = xB[0];
    const realtype xB1_ = xB[1];

    switch(ie) {
        case 0:
            switch(ip) {
                case 2:
                    deltaqB[0] = -xB0_;
                    break;
                case 3:
                    deltaqB[0] = xB1_;
                    break;
            }
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
