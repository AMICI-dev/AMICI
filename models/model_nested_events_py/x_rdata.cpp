#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"

namespace amici {
namespace model_model_nested_events_py {

void x_rdata_model_nested_events_py(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = Virus;
}

} // namespace model_model_nested_events_py
} // namespace amici
