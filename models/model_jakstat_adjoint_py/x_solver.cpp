#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_jakstat_adjoint_py {

void x_solver_model_jakstat_adjoint_py(realtype *x_solver, const realtype *x_rdata){
    const realtype STAT_ = x_rdata[0];
    const realtype pSTAT_ = x_rdata[1];
    const realtype pSTAT_pSTAT_ = x_rdata[2];
    const realtype npSTAT_npSTAT_ = x_rdata[3];
    const realtype nSTAT1_ = x_rdata[4];
    const realtype nSTAT2_ = x_rdata[5];
    const realtype nSTAT3_ = x_rdata[6];
    const realtype nSTAT4_ = x_rdata[7];
    const realtype nSTAT5_ = x_rdata[8];

    x_solver[0] = STAT_;
    x_solver[1] = pSTAT_;
    x_solver[2] = pSTAT_pSTAT_;
    x_solver[3] = npSTAT_npSTAT_;
    x_solver[4] = nSTAT1_;
    x_solver[5] = nSTAT2_;
    x_solver[6] = nSTAT3_;
    x_solver[7] = nSTAT4_;
    x_solver[8] = nSTAT5_;
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
