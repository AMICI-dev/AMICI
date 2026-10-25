#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

void rz_model_neuron_py(realtype *rz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h){
    const realtype v_ = x[0];

    switch(ie) {
        case 0:
            rz[0] = v_ - 30;
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
