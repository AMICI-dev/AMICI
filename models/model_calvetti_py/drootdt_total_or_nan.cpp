#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_calvetti_py {

void drootdt_total_or_nan_model_calvetti_py(realtype *drootdt_total_or_nan, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const realtype *tcl, const realtype *spl){
    drootdt_total_or_nan[0] = 1;
    drootdt_total_or_nan[1] = -1;
    drootdt_total_or_nan[2] = -1;
    drootdt_total_or_nan[3] = 1;
}

} // namespace model_model_calvetti_py
} // namespace amici
