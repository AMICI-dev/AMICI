#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_nested_events_py {

void x_rdata_model_nested_events_py(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    const realtype Virus_ = x[0];

    x_rdata[0] = Virus_;
}

} // namespace model_model_nested_events_py
} // namespace amici
