#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "p.h"
#include "k.h"
#include "rz.h"
#include "sigmaz.h"

namespace amici {
namespace model_model_neuron_py {

void dJrzdsigma_model_neuron_py(realtype *dJrzdsigma, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz){
    switch(iz) {
        case 0:
            dJrzdsigma[0] = -1.0*std::pow(rz1, 2)/std::pow(sigma_z1, 3) + 1.0/sigma_z1;
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
