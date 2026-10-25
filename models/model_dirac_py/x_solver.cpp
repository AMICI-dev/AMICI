#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_dirac_py {

void x_solver_model_dirac_py(realtype *x_solver, const realtype *x_rdata){
    const realtype x1_ = x_rdata[0];
    const realtype x2_ = x_rdata[1];

    x_solver[0] = x1_;
    x_solver[1] = x2_;
}

} // namespace model_model_dirac_py
} // namespace amici
