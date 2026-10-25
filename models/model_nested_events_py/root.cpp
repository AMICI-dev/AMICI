#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_nested_events_py {

void root_model_nested_events_py(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    const realtype Virus_ = x[0];
    const realtype t_0_ = p[2];

    root[0] = Virus_ - 1;
    root[1] = 1 - Virus_;
    root[2] = t - t_0_;
}

} // namespace model_model_nested_events_py
} // namespace amici
