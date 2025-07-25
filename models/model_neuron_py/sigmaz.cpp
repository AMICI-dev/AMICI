#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "sigmaz.h"

namespace amici {
namespace model_model_neuron_py {

void sigmaz_model_neuron_py(realtype *sigmaz, const realtype t, const realtype *p, const realtype *k){
    sigma_z1 = 1.0;  // sigmaz[0]
}

} // namespace model_model_neuron_py
} // namespace amici
