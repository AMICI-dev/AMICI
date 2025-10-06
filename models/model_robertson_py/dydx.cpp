#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_model_robertson_py {

void dydx_model_robertson_py(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dydx[0] = 1;
    dydx[4] = 10000.0;
    dydx[8] = 1;
}

} // namespace model_model_robertson_py
} // namespace amici
