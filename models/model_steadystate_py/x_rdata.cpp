#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_steadystate_py {

void x_rdata_model_steadystate_py(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    const realtype x1_ = x[0];
    const realtype x2_ = x[1];
    const realtype x3_ = x[2];

    x_rdata[0] = x1_;
    x_rdata[1] = x2_;
    x_rdata[2] = x3_;
}

} // namespace model_model_steadystate_py
} // namespace amici
