#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_dirac_py {

void deltaxB_at_trigger_model_dirac_py(realtype *deltaxB_at_trigger, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const int ie, const realtype *x_eval, const realtype *xB, const realtype *tcl){
    const realtype xB0_ = xB[0];

    switch(ie) {
        case 0:
            deltaxB_at_trigger[0] = xB0_;
            break;
    }
}

} // namespace model_model_dirac_py
} // namespace amici
