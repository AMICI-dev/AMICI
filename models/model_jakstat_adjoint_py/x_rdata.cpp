#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_jakstat_adjoint_py {

void x_rdata_model_jakstat_adjoint_py(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    const realtype STAT_ = x[0];
    const realtype pSTAT_ = x[1];
    const realtype pSTAT_pSTAT_ = x[2];
    const realtype npSTAT_npSTAT_ = x[3];
    const realtype nSTAT1_ = x[4];
    const realtype nSTAT2_ = x[5];
    const realtype nSTAT3_ = x[6];
    const realtype nSTAT4_ = x[7];
    const realtype nSTAT5_ = x[8];

    x_rdata[0] = STAT_;
    x_rdata[1] = pSTAT_;
    x_rdata[2] = pSTAT_pSTAT_;
    x_rdata[3] = npSTAT_npSTAT_;
    x_rdata[4] = nSTAT1_;
    x_rdata[5] = nSTAT2_;
    x_rdata[6] = nSTAT3_;
    x_rdata[7] = nSTAT4_;
    x_rdata[8] = nSTAT5_;
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
