#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x_rdata.h"

namespace amici {
namespace model_model_neuron_py {

void x_solver_model_neuron_py(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = v;
    x_solver[1] = u;
}

} // namespace model_model_neuron_py
} // namespace amici
