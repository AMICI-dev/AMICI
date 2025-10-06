#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "k.h"
#include "y.h"
#include "sigmay.h"
#include "my.h"

namespace amici {
namespace model_model_calvetti_py {

void dJydsigma_model_calvetti_py(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_obs_V1 - 1.0*std::pow(-mobs_V1 + obs_V1, 2)/std::pow(sigma_obs_V1, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_obs_V2 - 1.0*std::pow(-mobs_V2 + obs_V2, 2)/std::pow(sigma_obs_V2, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_obs_V3 - 1.0*std::pow(-mobs_V3 + obs_V3, 2)/std::pow(sigma_obs_V3, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_obs_f0 - 1.0*std::pow(-mobs_f0 + obs_f0, 2)/std::pow(sigma_obs_f0, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_obs_f1 - 1.0*std::pow(-mobs_f1 + obs_f1, 2)/std::pow(sigma_obs_f1, 3);
            break;
        case 5:
            dJydsigma[5] = 1.0/sigma_obs_f2 - 1.0*std::pow(-mobs_f2 + obs_f2, 2)/std::pow(sigma_obs_f2, 3);
            break;
    }
}

} // namespace model_model_calvetti_py
} // namespace amici
