#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_nested_events_py {

void deltax_model_nested_events_py(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_eval, const realtype *x_eval){
    const realtype Virus_ = x[0];
    const realtype V_0_inject_ = p[1];
    const realtype x_eval0_ = x_eval[0];

    switch(ie) {
        case 2:
            deltax[0] = V_0_inject_ - Virus_ + x_eval0_;
            break;
    }
}

} // namespace model_model_nested_events_py
} // namespace amici
