#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_events_py {

void dydp_model_events_py(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl){
    const realtype x1_ = x[0];
    const realtype x2_ = x[1];
    const realtype x3_ = x[2];

    switch(ip) {
        case 3:
            dydp[0] = x1_ + x2_ + x3_;
            break;
    }
}

} // namespace model_model_events_py
} // namespace amici
