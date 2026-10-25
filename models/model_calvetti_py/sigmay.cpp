#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_calvetti_py {

void sigmay_model_calvetti_py(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    realtype &sigma_obs_V1_ = sigmay[0];
    realtype &sigma_obs_V2_ = sigmay[1];
    realtype &sigma_obs_V3_ = sigmay[2];
    realtype &sigma_obs_f0_ = sigmay[3];
    realtype &sigma_obs_f1_ = sigmay[4];
    realtype &sigma_obs_f2_ = sigmay[5];
    sigma_obs_V1_ = 1.0;  // sigmay[0]
    sigma_obs_V2_ = 1.0;  // sigmay[1]
    sigma_obs_V3_ = 1.0;  // sigmay[2]
    sigma_obs_f0_ = 1.0;  // sigmay[3]
    sigma_obs_f1_ = 1.0;  // sigmay[4]
    sigma_obs_f2_ = 1.0;  // sigmay[5]
}

} // namespace model_model_calvetti_py
} // namespace amici
