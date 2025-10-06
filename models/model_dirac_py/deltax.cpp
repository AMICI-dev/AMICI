#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "h.h"
#include "xdot.h"
#include "x_old.h"

namespace amici {
namespace model_model_dirac_py {

void deltax_model_dirac_py(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old){
    switch(ie) {
        case 0:
            deltax[0] = -x1 + x_old0 + 1;
            break;
    }
}

} // namespace model_model_dirac_py
} // namespace amici
