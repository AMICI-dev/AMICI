#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_events_py {

void drootdt_total_or_nan_model_events_py(realtype *drootdt_total_or_nan, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const realtype *tcl, const realtype *spl){
    const realtype x1_ = x[0];
    const realtype x2_ = x[1];
    const realtype x3_ = x[2];
    const realtype p1_ = p[0];
    const realtype p2_ = p[1];
    const realtype p3_ = p[2];
    const realtype Heaviside_2_ = h[2];
    const realtype Heaviside_4_ = h[4];

    drootdt_total_or_nan[0] = Heaviside_4_ + p2_*x1_*std::exp(-1.0/10.0*t) - p3_*x2_ + x3_ - 1;
    drootdt_total_or_nan[1] = Heaviside_4_ - p1_*x1_*(1 - Heaviside_2_) + x3_ - 1;
    drootdt_total_or_nan[2] = -1;
    drootdt_total_or_nan[3] = 1;
    drootdt_total_or_nan[4] = -1;
    drootdt_total_or_nan[5] = 1;
}

} // namespace model_model_events_py
} // namespace amici
