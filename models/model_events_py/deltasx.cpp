#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"
#include "xdot.h"
#include "xdot_old.h"
#include "sx.h"
#include "stau.h"
#include "x_old.h"

namespace amici {
namespace model_model_events_py {

void deltasx_model_events_py(realtype *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau, const realtype *tcl, const realtype *x_old){
    switch(ie) {
        case 0:
        case 1:
        case 2:
        case 3:
            switch(ip) {
                case 0:
                case 1:
                case 2:
                case 3:
                    deltasx[0] = stau0*(dx1dt - xdot_old0);
                    deltasx[1] = stau0*(dx2dt - xdot_old1);
                    deltasx[2] = stau0*(dx3dt - xdot_old2);
                    break;
            }
            break;
    }
}

} // namespace model_model_events_py
} // namespace amici
