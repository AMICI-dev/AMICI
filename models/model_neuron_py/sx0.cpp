#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

void sx0_model_neuron_py(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip){
    const realtype v0_ = k[0];

    switch(ip) {
        case 1:
            sx0[1] = v0_;
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
