#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

void deltaxB_model_neuron_py(realtype *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const int ie, const realtype *x_eval, const realtype *xB, const realtype *tcl){
    const realtype xB0_ = xB[0];

    switch(ie) {
        case 0:
            deltaxB[0] = -xB0_;
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
