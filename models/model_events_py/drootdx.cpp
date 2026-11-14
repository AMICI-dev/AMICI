#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_events_py {

void drootdx_model_events_py(realtype *drootdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const realtype *tcl, const realtype *spl){
    drootdx[1] = 1;
    drootdx[2] = -1;
    drootdx[3] = 1;
    drootdx[5] = -1;
}

} // namespace model_model_events_py
} // namespace amici
