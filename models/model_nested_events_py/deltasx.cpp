#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_nested_events_py {

void deltasx_model_nested_events_py(realtype *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_eval, const realtype *sx, const realtype *stau, const realtype *tcl, const realtype *x_eval){
    const realtype dVirusdt_ = xdot[0];
    const realtype xdot_eval0_ = xdot_eval[0];
    const realtype stau0_ = stau[0];

    switch(ie) {
        case 0:
        case 1:
            switch(ip) {
                case 0:
                case 1:
                case 2:
                case 3:
                case 4:
                    deltasx[0] = stau0_*(dVirusdt_ - xdot_eval0_);
                    break;
            }
            break;
        case 2:
            switch(ip) {
                case 0:
                case 2:
                case 3:
                case 4:
                    deltasx[0] = stau0_*(dVirusdt_ - xdot_eval0_);
                    break;
                case 1:
                    deltasx[0] = stau0_*(dVirusdt_ - xdot_eval0_) + 1;
                    break;
            }
            break;
    }
}

} // namespace model_model_nested_events_py
} // namespace amici
