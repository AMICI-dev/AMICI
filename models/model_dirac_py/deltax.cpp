#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_dirac_py {

void deltax_model_dirac_py(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_eval, const realtype *x_eval){
    const realtype x1_ = x[0];
    const realtype x_eval0_ = x_eval[0];

    switch(ie) {
        case 0:
            deltax[0] = -x1_ + x_eval0_ + 1;
            break;
    }
}

} // namespace model_model_dirac_py
} // namespace amici
