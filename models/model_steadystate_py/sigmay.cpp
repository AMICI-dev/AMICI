#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_steadystate_py {

void sigmay_model_steadystate_py(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    realtype &sigma_obs_x1_ = sigmay[0];
    realtype &sigma_obs_x2_ = sigmay[1];
    realtype &sigma_obs_x3_ = sigmay[2];
    sigma_obs_x1_ = 1.0;  // sigmay[0]
    sigma_obs_x2_ = 1.0;  // sigmay[1]
    sigma_obs_x3_ = 1.0;  // sigmay[2]
}

} // namespace model_model_steadystate_py
} // namespace amici
