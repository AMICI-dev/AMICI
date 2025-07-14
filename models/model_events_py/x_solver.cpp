#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_model_events_py {

void x_solver_model_events_py(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = x1;
    x_solver[1] = x2;
    x_solver[2] = x3;
}

} // namespace model_model_events_py
} // namespace amici
