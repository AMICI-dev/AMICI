#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "k.h"

namespace amici {
namespace model_model_calvetti_py {

void x_rdata_model_calvetti_py(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = V1;
    x_rdata[1] = V2;
    x_rdata[2] = V3;
    x_rdata[3] = f1;
    x_rdata[4] = f2;
    x_rdata[5] = f3;
}

} // namespace model_model_calvetti_py
} // namespace amici
