#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_nested_events_py {

void deltaqB_model_nested_events_py(realtype *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const int ip, const int ie, const realtype *x_eval, const realtype *xB){
    const realtype xB0_ = xB[0];

    switch(ie) {
        case 2:
            switch(ip) {
                case 1:
                    deltaqB[0] = xB0_;
                    break;
            }
            break;
    }
}

} // namespace model_model_nested_events_py
} // namespace amici
