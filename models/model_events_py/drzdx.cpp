#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"

namespace amici {
namespace model_model_events_py {

void drzdx_model_events_py(realtype *drzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h){
    switch(ie) {
        case 0:
            drzdx[2] = 1;
            drzdx[4] = -1;
            break;
        case 1:
            drzdx[1] = 1;
            drzdx[5] = -1;
            break;
    }
}

} // namespace model_model_events_py
} // namespace amici
