#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_dirac_py {

void x_rdata_model_dirac_py(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    const realtype x1_ = x[0];
    const realtype x2_ = x[1];

    x_rdata[0] = x1_;
    x_rdata[1] = x2_;
}

} // namespace model_model_dirac_py
} // namespace amici
