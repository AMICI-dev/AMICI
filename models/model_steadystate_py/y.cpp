#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_steadystate_py {

void y_model_steadystate_py(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype x1_ = x[0];
    const realtype x2_ = x[1];
    const realtype x3_ = x[2];

    y[0] = x1_;
    y[1] = x2_;
    y[2] = x3_;
}

} // namespace model_model_steadystate_py
} // namespace amici
