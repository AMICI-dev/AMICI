#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

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
