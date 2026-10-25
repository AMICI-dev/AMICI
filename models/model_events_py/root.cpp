#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_events_py {

void root_model_events_py(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    const realtype x1_ = x[0];
    const realtype x2_ = x[1];
    const realtype x3_ = x[2];
    const realtype p4_ = p[3];

    root[0] = x2_ - x3_;
    root[1] = x1_ - x3_;
    root[2] = p4_ - t;
    root[3] = -p4_ + t;
    root[4] = 4 - t;
    root[5] = t - 4;
}

} // namespace model_model_events_py
} // namespace amici
