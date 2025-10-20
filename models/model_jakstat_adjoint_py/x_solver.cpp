#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x_rdata.h"

namespace amici {
namespace model_model_jakstat_adjoint_py {

void x_solver_model_jakstat_adjoint_py(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = STAT;
    x_solver[1] = pSTAT;
    x_solver[2] = pSTAT_pSTAT;
    x_solver[3] = npSTAT_npSTAT;
    x_solver[4] = nSTAT1;
    x_solver[5] = nSTAT2;
    x_solver[6] = nSTAT3;
    x_solver[7] = nSTAT4;
    x_solver[8] = nSTAT5;
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
