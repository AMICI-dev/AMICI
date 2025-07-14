#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_model_jakstat_adjoint_py {

void x_rdata_model_jakstat_adjoint_py(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = STAT;
    x_rdata[1] = pSTAT;
    x_rdata[2] = pSTAT_pSTAT;
    x_rdata[3] = npSTAT_npSTAT;
    x_rdata[4] = nSTAT1;
    x_rdata[5] = nSTAT2;
    x_rdata[6] = nSTAT3;
    x_rdata[7] = nSTAT4;
    x_rdata[8] = nSTAT5;
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
