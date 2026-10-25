#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_nested_events_py {

void x_solver_model_nested_events_py(realtype *x_solver, const realtype *x_rdata){
    const realtype Virus_ = x_rdata[0];

    x_solver[0] = Virus_;
}

} // namespace model_model_nested_events_py
} // namespace amici
