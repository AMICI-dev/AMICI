#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_dirac_py {

void drootdt_total_or_nan_model_dirac_py(realtype *drootdt_total_or_nan, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const realtype *tcl, const realtype *spl){
    drootdt_total_or_nan[0] = 1;
}

} // namespace model_model_dirac_py
} // namespace amici
