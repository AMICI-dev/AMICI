#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "h.h"

namespace amici {
namespace model_model_dirac_py {

void root_model_dirac_py(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    root[0] = -p2 + t;
}

} // namespace model_model_dirac_py
} // namespace amici
