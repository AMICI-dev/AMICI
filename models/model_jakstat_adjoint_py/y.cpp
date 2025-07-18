#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"

namespace amici {
namespace model_model_jakstat_adjoint_py {

void y_model_jakstat_adjoint_py(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = offset_pSTAT + scale_pSTAT*(pSTAT + 2*pSTAT_pSTAT)/init_STAT;
    y[1] = offset_tSTAT + scale_tSTAT*(STAT + pSTAT + 2*pSTAT_pSTAT)/init_STAT;
    y[2] = u;
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
