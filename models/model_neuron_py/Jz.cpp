#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

void Jz_model_neuron_py(realtype *Jz, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const realtype *mz){
    const realtype z1_ = z[0];
    const realtype sigma_z1_ = sigmaz[0];
    const realtype mz1_ = mz[0];

    switch(iz) {
        case 0:
            Jz[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_z1_, 2)) + 0.5*std::pow(-mz1_ + z1_, 2)/std::pow(sigma_z1_, 2);
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
