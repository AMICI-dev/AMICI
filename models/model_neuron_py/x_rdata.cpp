#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

void x_rdata_model_neuron_py(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    const realtype v_ = x[0];
    const realtype u_ = x[1];

    x_rdata[0] = v_;
    x_rdata[1] = u_;
}

} // namespace model_model_neuron_py
} // namespace amici
