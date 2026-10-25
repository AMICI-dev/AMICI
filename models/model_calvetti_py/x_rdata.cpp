#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_calvetti_py {

void x_rdata_model_calvetti_py(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    const realtype V1_ = x[0];
    const realtype V2_ = x[1];
    const realtype V3_ = x[2];
    const realtype f1_ = x[3];
    const realtype f2_ = x[4];
    const realtype f3_ = x[5];

    x_rdata[0] = V1_;
    x_rdata[1] = V2_;
    x_rdata[2] = V3_;
    x_rdata[3] = f1_;
    x_rdata[4] = f2_;
    x_rdata[5] = f3_;
}

} // namespace model_model_calvetti_py
} // namespace amici
