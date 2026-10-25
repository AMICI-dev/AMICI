#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_jakstat_adjoint_py {

void sigmay_model_jakstat_adjoint_py(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    const realtype sigma_pSTAT_ = p[14];
    const realtype sigma_tSTAT_ = p[15];
    const realtype sigma_pEpoR_ = p[16];

    realtype &sigma_obs_pSTAT_ = sigmay[0];
    realtype &sigma_obs_tSTAT_ = sigmay[1];
    realtype &sigma_obs_spline_ = sigmay[2];
    sigma_obs_pSTAT_ = sigma_pSTAT_;  // sigmay[0]
    sigma_obs_tSTAT_ = sigma_tSTAT_;  // sigmay[1]
    sigma_obs_spline_ = sigma_pEpoR_;  // sigmay[2]
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
