#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_robertson_py {

void x_solver_model_robertson_py(realtype *x_solver, const realtype *x_rdata){
    const realtype x1_ = x_rdata[0];
    const realtype x2_ = x_rdata[1];
    const realtype x3_ = x_rdata[2];

    x_solver[0] = x1_;
    x_solver[1] = x2_;
    x_solver[2] = x3_;
}

} // namespace model_model_robertson_py
} // namespace amici
