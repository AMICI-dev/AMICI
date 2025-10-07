#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"
#include "sx.h"

namespace amici {
namespace model_model_events_py {

void stau_model_events_py(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *tcl, const realtype *sx, const int ip, const int ie){
    switch(ie) {
        case 0:
            switch(ip) {
                case 0:
                case 1:
                case 2:
                case 3:
                    stau[0] = (sx1 - sx2)/(Heaviside_4 + p2*x1*std::exp(-1.0/10.0*t) - p3*x2 + x3 - 1);
                    break;
            }
            break;
        case 1:
            switch(ip) {
                case 0:
                case 1:
                case 2:
                case 3:
                    stau[0] = (sx0 - sx2)/(Heaviside_4 - p1*x1*(1 - Heaviside_2) + x3 - 1);
                    break;
            }
            break;
        case 2:
        case 3:
            switch(ip) {
                case 3:
                    stau[0] = -1;
                    break;
            }
            break;
    }
}

} // namespace model_model_events_py
} // namespace amici
