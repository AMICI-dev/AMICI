#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"

namespace amici {
namespace model_model_dirac_py {

void x_rdata_model_dirac_py(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = x1;
    x_rdata[1] = x2;
}

} // namespace model_model_dirac_py
} // namespace amici
