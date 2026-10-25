#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_events_py {

void y_model_events_py(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype x1_ = x[0];
    const realtype x2_ = x[1];
    const realtype x3_ = x[2];
    const realtype p4_ = p[3];

    y[0] = p4_*(x1_ + x2_ + x3_);
}

} // namespace model_model_events_py
} // namespace amici
