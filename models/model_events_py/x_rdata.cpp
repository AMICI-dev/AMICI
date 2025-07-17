#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_model_events_py {

void x_rdata_model_events_py(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = x1;
    x_rdata[1] = x2;
    x_rdata[2] = x3;
}

} // namespace model_model_events_py
} // namespace amici
