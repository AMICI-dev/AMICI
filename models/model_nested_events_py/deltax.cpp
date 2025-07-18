#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "h.h"
#include "xdot.h"
#include "x_old.h"

namespace amici {
namespace model_model_nested_events_py {

void deltax_model_nested_events_py(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old){
    switch(ie) {
        case 2:
            deltax[0] = V_0_inject - Virus + x_old0;
            break;
    }
}

} // namespace model_model_nested_events_py
} // namespace amici
