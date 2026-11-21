#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_nested_events_py {

void drootdt_total_or_nan_model_nested_events_py(realtype *drootdt_total_or_nan, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const realtype *tcl, const realtype *spl){
    const realtype Virus_ = x[0];
    const realtype rho_V_ = p[3];
    const realtype delta_V_ = p[4];
    const realtype Heaviside_1_ = h[0];

    drootdt_total_or_nan[0] = Heaviside_1_*Virus_*rho_V_ - Virus_*delta_V_;
    drootdt_total_or_nan[1] = -Heaviside_1_*Virus_*rho_V_ + Virus_*delta_V_;
    drootdt_total_or_nan[2] = 1;
}

} // namespace model_model_nested_events_py
} // namespace amici
