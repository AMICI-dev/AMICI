#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_events_py {

void stau_model_events_py(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const realtype *tcl, const realtype *sx, const int ip, const int ie){
    const realtype x1_ = x[0];
    const realtype x2_ = x[1];
    const realtype x3_ = x[2];
    const realtype p1_ = p[0];
    const realtype p2_ = p[1];
    const realtype p3_ = p[2];
    const realtype Heaviside_2_ = h[2];
    const realtype Heaviside_4_ = h[4];
    const realtype sx0_ = sx[0];
    const realtype sx1_ = sx[1];
    const realtype sx2_ = sx[2];

    switch(ie) {
        case 0:
            switch(ip) {
                case 0:
                case 1:
                case 2:
                case 3:
                    stau[0] = (sx1_ - sx2_)/(Heaviside_4_ + p2_*x1_*std::exp(-1.0/10.0*t) - p3_*x2_ + x3_ - 1);
                    break;
            }
            break;
        case 1:
            switch(ip) {
                case 0:
                case 1:
                case 2:
                case 3:
                    stau[0] = (sx0_ - sx2_)/(Heaviside_4_ - p1_*x1_*(1 - Heaviside_2_) + x3_ - 1);
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
