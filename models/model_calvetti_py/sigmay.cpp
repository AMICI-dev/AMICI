#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "k.h"
#include "y.h"
#include "sigmay.h"

namespace amici {
namespace model_model_calvetti_py {

void sigmay_model_calvetti_py(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_obs_V1 = 1.0;  // sigmay[0]
    sigma_obs_V2 = 1.0;  // sigmay[1]
    sigma_obs_V3 = 1.0;  // sigmay[2]
    sigma_obs_f0 = 1.0;  // sigmay[3]
    sigma_obs_f1 = 1.0;  // sigmay[4]
    sigma_obs_f2 = 1.0;  // sigmay[5]
}

} // namespace model_model_calvetti_py
} // namespace amici
