#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "k.h"
#include "h.h"
#include "w.h"

namespace amici {
namespace model_model_calvetti_py {

void y_model_calvetti_py(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = V1;
    y[1] = V2;
    y[2] = V3;
    y[3] = f0;
    y[4] = f1;
    y[5] = f2;
}

} // namespace model_model_calvetti_py
} // namespace amici
