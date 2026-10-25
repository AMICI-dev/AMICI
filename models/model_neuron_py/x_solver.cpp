#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

void x_solver_model_neuron_py(realtype *x_solver, const realtype *x_rdata){
    const realtype v_ = x_rdata[0];
    const realtype u_ = x_rdata[1];

    x_solver[0] = v_;
    x_solver[1] = u_;
}

} // namespace model_model_neuron_py
} // namespace amici
