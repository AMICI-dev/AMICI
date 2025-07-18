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

void y_model_events_py(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = p4*(x1 + x2 + x3);
}

} // namespace model_model_events_py
} // namespace amici
