#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

void dJrzdrz_model_neuron_py(realtype *dJrzdrz, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz){
    const realtype rz1_ = rz[0];
    const realtype sigma_z1_ = sigmaz[0];

    switch(iz) {
        case 0:
            dJrzdrz[0] = 1.0*rz1_/std::pow(sigma_z1_, 2);
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
