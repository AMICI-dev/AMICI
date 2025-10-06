#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x_rdata.h"

namespace amici {
namespace model_model_calvetti_py {

void x_solver_model_calvetti_py(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = V1;
    x_solver[1] = V2;
    x_solver[2] = V3;
    x_solver[3] = f1;
    x_solver[4] = f2;
    x_solver[5] = f3;
}

} // namespace model_model_calvetti_py
} // namespace amici
