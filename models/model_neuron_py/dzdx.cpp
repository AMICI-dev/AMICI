#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

void dzdx_model_neuron_py(realtype *dzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h){
    const realtype v_ = x[0];
    const realtype u_ = x[1];
    const realtype I0_ = k[1];

    switch(ie) {
        case 0:
            dzdx[0] = -1/(I0_ - u_ + (1.0/25.0)*std::pow(v_, 2) + 5*v_ + 140);
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
