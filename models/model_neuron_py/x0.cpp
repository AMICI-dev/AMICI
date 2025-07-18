#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_model_neuron_py {

void x0_model_neuron_py(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = v0;
    x0[1] = b*v0;
}

} // namespace model_model_neuron_py
} // namespace amici
