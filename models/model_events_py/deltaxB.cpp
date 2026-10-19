#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_events_py {

void deltaxB_model_events_py(realtype *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const int ie, const realtype *xdot, const realtype *xdot_eval, const realtype *x_eval, const realtype *xB, const realtype *tcl){
    const realtype p1_ = p[0];
    const realtype p2_ = p[1];
    const realtype p3_ = p[2];
    const realtype Heaviside_2_ = h[2];
    const realtype Heaviside_4_ = h[4];
    const realtype dx1dt_ = xdot[0];
    const realtype dx2dt_ = xdot[1];
    const realtype dx3dt_ = xdot[2];
    const realtype xdot_eval0_ = xdot_eval[0];
    const realtype xdot_eval1_ = xdot_eval[1];
    const realtype xdot_eval2_ = xdot_eval[2];
    const realtype x_eval0_ = x_eval[0];
    const realtype x_eval1_ = x_eval[1];
    const realtype x_eval2_ = x_eval[2];
    const realtype xB0_ = xB[0];
    const realtype xB1_ = xB[1];
    const realtype xB2_ = xB[2];

    switch(ie) {
        case 0:
            deltaxB[1] = xB0_*(dx1dt_ - xdot_eval0_)/(Heaviside_4_ + p2_*x_eval0_*std::exp(-1.0/10.0*t) - p3_*x_eval1_ + x_eval2_ - 1) + xB1_*(dx2dt_ - xdot_eval1_)/(Heaviside_4_ + p2_*x_eval0_*std::exp(-1.0/10.0*t) - p3_*x_eval1_ + x_eval2_ - 1) + xB2_*(dx3dt_ - xdot_eval2_)/(Heaviside_4_ + p2_*x_eval0_*std::exp(-1.0/10.0*t) - p3_*x_eval1_ + x_eval2_ - 1);
            deltaxB[2] = -xB0_*(dx1dt_ - xdot_eval0_)/(Heaviside_4_ + p2_*x_eval0_*std::exp(-1.0/10.0*t) - p3_*x_eval1_ + x_eval2_ - 1) - xB1_*(dx2dt_ - xdot_eval1_)/(Heaviside_4_ + p2_*x_eval0_*std::exp(-1.0/10.0*t) - p3_*x_eval1_ + x_eval2_ - 1) - xB2_*(dx3dt_ - xdot_eval2_)/(Heaviside_4_ + p2_*x_eval0_*std::exp(-1.0/10.0*t) - p3_*x_eval1_ + x_eval2_ - 1);
            break;
        case 1:
            deltaxB[0] = xB0_*(dx1dt_ - xdot_eval0_)/(Heaviside_4_ - p1_*x_eval0_*(1 - Heaviside_2_) + x_eval2_ - 1) + xB1_*(dx2dt_ - xdot_eval1_)/(Heaviside_4_ - p1_*x_eval0_*(1 - Heaviside_2_) + x_eval2_ - 1) + xB2_*(dx3dt_ - xdot_eval2_)/(Heaviside_4_ - p1_*x_eval0_*(1 - Heaviside_2_) + x_eval2_ - 1);
            deltaxB[2] = -xB0_*(dx1dt_ - xdot_eval0_)/(Heaviside_4_ - p1_*x_eval0_*(1 - Heaviside_2_) + x_eval2_ - 1) - xB1_*(dx2dt_ - xdot_eval1_)/(Heaviside_4_ - p1_*x_eval0_*(1 - Heaviside_2_) + x_eval2_ - 1) - xB2_*(dx3dt_ - xdot_eval2_)/(Heaviside_4_ - p1_*x_eval0_*(1 - Heaviside_2_) + x_eval2_ - 1);
            break;
    }
}

} // namespace model_model_events_py
} // namespace amici
