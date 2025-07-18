#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_model_robertson_py {

void y_model_robertson_py(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = x1;
    y[1] = 10000.0*x2;
    y[2] = x3;
}

} // namespace model_model_robertson_py
} // namespace amici
