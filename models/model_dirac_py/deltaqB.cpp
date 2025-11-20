#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "h.h"
#include "xdot.h"
#include "xdot_old.h"
#include "x_old.h"
#include "xB.h"

namespace amici {
namespace model_model_dirac_py {

void deltaqB_model_dirac_py(realtype *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old, const realtype *xB){
    switch(ie) {
        case 0:
            switch(ip) {
                case 1:
                    deltaqB[0] = xB0*(-dx1dt + xdot_old0) + xB1*(-dx2dt + xdot_old1);
                    break;
            }
            break;
    }
}

} // namespace model_model_dirac_py
} // namespace amici
