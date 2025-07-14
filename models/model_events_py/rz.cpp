#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"
#include "rz.h"

namespace amici {
namespace model_model_events_py {

void rz_model_events_py(realtype *rz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h){
    switch(ie) {
        case 0:
            rz[0] = x2 - x3;
            break;
        case 1:
            rz[1] = x1 - x3;
            break;
    }
}

} // namespace model_model_events_py
} // namespace amici
