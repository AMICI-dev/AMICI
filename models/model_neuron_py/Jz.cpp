#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "p.h"
#include "k.h"
#include "z.h"
#include "sigmaz.h"
#include "mz.h"

namespace amici {
namespace model_model_neuron_py {

void Jz_model_neuron_py(realtype *Jz, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const realtype *mz){
    switch(iz) {
        case 0:
            Jz[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_z1, 2)) + 0.5*std::pow(-mz1 + z1, 2)/std::pow(sigma_z1, 2);
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
