#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "k.h"
#include "y.h"
#include "sigmay.h"
#include "my.h"

namespace amici {
namespace model_model_calvetti_py {

void Jy_model_calvetti_py(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_V1, 2)) + 0.5*std::pow(-mobs_V1 + obs_V1, 2)/std::pow(sigma_obs_V1, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_V2, 2)) + 0.5*std::pow(-mobs_V2 + obs_V2, 2)/std::pow(sigma_obs_V2, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_V3, 2)) + 0.5*std::pow(-mobs_V3 + obs_V3, 2)/std::pow(sigma_obs_V3, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_f0, 2)) + 0.5*std::pow(-mobs_f0 + obs_f0, 2)/std::pow(sigma_obs_f0, 2);
            break;
        case 4:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_f1, 2)) + 0.5*std::pow(-mobs_f1 + obs_f1, 2)/std::pow(sigma_obs_f1, 2);
            break;
        case 5:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_f2, 2)) + 0.5*std::pow(-mobs_f2 + obs_f2, 2)/std::pow(sigma_obs_f2, 2);
            break;
    }
}

} // namespace model_model_calvetti_py
} // namespace amici
