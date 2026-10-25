#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_dirac_py {

void y_model_dirac_py(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype x2_ = x[1];

    y[0] = x2_;
}

} // namespace model_model_dirac_py
} // namespace amici
