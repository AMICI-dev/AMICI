#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "spl.h"
#include "w.h"

namespace amici {
namespace model_model_jakstat_adjoint_py {

void w_model_jakstat_adjoint_py(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){

    // dynamic expressions
    u = spl_0;  // w[0]
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
