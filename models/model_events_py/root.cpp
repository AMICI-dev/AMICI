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

void root_model_events_py(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    root[0] = x2 - x3;
    root[1] = x1 - x3;
    root[2] = p4 - t;
    root[3] = -p4 + t;
    root[4] = 4 - t;
    root[5] = t - 4;
}

} // namespace model_model_events_py
} // namespace amici
