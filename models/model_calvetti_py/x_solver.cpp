#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_calvetti_py {

void x_solver_model_calvetti_py(realtype *x_solver, const realtype *x_rdata){
    const realtype V1_ = x_rdata[0];
    const realtype V2_ = x_rdata[1];
    const realtype V3_ = x_rdata[2];
    const realtype f1_ = x_rdata[3];
    const realtype f2_ = x_rdata[4];
    const realtype f3_ = x_rdata[5];

    x_solver[0] = V1_;
    x_solver[1] = V2_;
    x_solver[2] = V3_;
    x_solver[3] = f1_;
    x_solver[4] = f2_;
    x_solver[5] = f3_;
}

} // namespace model_model_calvetti_py
} // namespace amici
