#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "h.h"
#include "xdot.h"
#include "xdot_old.h"
#include "sx.h"
#include "stau.h"
#include "x_old.h"

namespace amici {
namespace model_model_nested_events_py {

void deltasx_model_nested_events_py(realtype *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau, const realtype *tcl, const realtype *x_old){
    switch(ie) {
        case 0:
        case 1:
            switch(ip) {
                case 0:
                case 1:
                case 2:
                case 3:
                case 4:
                    deltasx[0] = stau0*(dVirusdt - xdot_old0);
                    break;
            }
            break;
        case 2:
            switch(ip) {
                case 0:
                case 2:
                case 3:
                case 4:
                    deltasx[0] = stau0*(dVirusdt - xdot_old0);
                    break;
                case 1:
                    deltasx[0] = stau0*(dVirusdt - xdot_old0) + 1;
                    break;
            }
            break;
    }
}

} // namespace model_model_nested_events_py
} // namespace amici
