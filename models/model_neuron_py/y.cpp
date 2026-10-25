#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

void y_model_neuron_py(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype v_ = x[0];

    y[0] = v_;
}

} // namespace model_model_neuron_py
} // namespace amici
