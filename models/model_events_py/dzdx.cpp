#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"

namespace amici {
namespace model_model_events_py {

void dzdx_model_events_py(realtype *dzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h){
    switch(ie) {
        case 0:
            dzdx[2] = -1/(Heaviside_4 + p2*x1*std::exp(-1.0/10.0*t) - p3*x2 + x3 - 1);
            dzdx[4] = 1.0/(Heaviside_4 + p2*x1*std::exp(-1.0/10.0*t) - p3*x2 + x3 - 1);
            break;
        case 1:
            dzdx[1] = -1/(Heaviside_4 - p1*x1*(1 - Heaviside_2) + x3 - 1);
            dzdx[5] = 1.0/(Heaviside_4 - p1*x1*(1 - Heaviside_2) + x3 - 1);
            break;
    }
}

} // namespace model_model_events_py
} // namespace amici
