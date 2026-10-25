#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

void sigmaz_model_neuron_py(realtype *sigmaz, const realtype t, const realtype *p, const realtype *k){
    realtype &sigma_z1_ = sigmaz[0];
    sigma_z1_ = 1.0;  // sigmaz[0]
}

} // namespace model_model_neuron_py
} // namespace amici
