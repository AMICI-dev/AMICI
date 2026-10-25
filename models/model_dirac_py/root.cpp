#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_dirac_py {

void root_model_dirac_py(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    const realtype p2_ = p[1];

    root[0] = -p2_ + t;
}

} // namespace model_model_dirac_py
} // namespace amici
