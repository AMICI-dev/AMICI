#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

void root_model_neuron_py(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    const realtype v_ = x[0];

    root[0] = v_ - 30;
}

} // namespace model_model_neuron_py
} // namespace amici
