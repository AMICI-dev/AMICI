#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"
#include "sigmay.h"

namespace amici {
namespace model_model_jakstat_adjoint_py {

void sigmay_model_jakstat_adjoint_py(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_obs_pSTAT = sigma_pSTAT;  // sigmay[0]
    sigma_obs_tSTAT = sigma_tSTAT;  // sigmay[1]
    sigma_obs_spline = sigma_pEpoR;  // sigmay[2]
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
