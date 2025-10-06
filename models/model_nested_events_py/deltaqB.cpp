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
namespace model_model_nested_events_py {

void deltaqB_model_nested_events_py(realtype *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *x_old, const realtype *xB){
    switch(ie) {
        case 2:
            switch(ip) {
                case 1:
                    deltaqB[0] = xB0;
                    break;
                case 2:
                    deltaqB[0] = xB0*(-dVirusdt + xdot_old0);
                    break;
            }
            break;
    }
}

} // namespace model_model_nested_events_py
} // namespace amici
