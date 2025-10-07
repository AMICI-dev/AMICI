#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "h.h"

namespace amici {
namespace model_model_nested_events_py {

void root_model_nested_events_py(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    root[0] = Virus - 1;
    root[1] = 1 - Virus;
    root[2] = t - t_0;
}

} // namespace model_model_nested_events_py
} // namespace amici
