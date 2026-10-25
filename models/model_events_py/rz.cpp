#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_events_py {

void rz_model_events_py(realtype *rz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h){
    const realtype x1_ = x[0];
    const realtype x2_ = x[1];
    const realtype x3_ = x[2];

    switch(ie) {
        case 0:
            rz[0] = x2_ - x3_;
            break;
        case 1:
            rz[1] = x1_ - x3_;
            break;
    }
}

} // namespace model_model_events_py
} // namespace amici
