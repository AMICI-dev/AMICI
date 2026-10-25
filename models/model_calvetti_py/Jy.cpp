#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_calvetti_py {

void Jy_model_calvetti_py(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    const realtype obs_V1_ = y[0];
    const realtype obs_V2_ = y[1];
    const realtype obs_V3_ = y[2];
    const realtype obs_f0_ = y[3];
    const realtype obs_f1_ = y[4];
    const realtype obs_f2_ = y[5];
    const realtype sigma_obs_V1_ = sigmay[0];
    const realtype sigma_obs_V2_ = sigmay[1];
    const realtype sigma_obs_V3_ = sigmay[2];
    const realtype sigma_obs_f0_ = sigmay[3];
    const realtype sigma_obs_f1_ = sigmay[4];
    const realtype sigma_obs_f2_ = sigmay[5];
    const realtype mobs_V1_ = my[0];
    const realtype mobs_V2_ = my[1];
    const realtype mobs_V3_ = my[2];
    const realtype mobs_f0_ = my[3];
    const realtype mobs_f1_ = my[4];
    const realtype mobs_f2_ = my[5];

    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_V1_, 2)) + 0.5*std::pow(-mobs_V1_ + obs_V1_, 2)/std::pow(sigma_obs_V1_, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_V2_, 2)) + 0.5*std::pow(-mobs_V2_ + obs_V2_, 2)/std::pow(sigma_obs_V2_, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_V3_, 2)) + 0.5*std::pow(-mobs_V3_ + obs_V3_, 2)/std::pow(sigma_obs_V3_, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_f0_, 2)) + 0.5*std::pow(-mobs_f0_ + obs_f0_, 2)/std::pow(sigma_obs_f0_, 2);
            break;
        case 4:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_f1_, 2)) + 0.5*std::pow(-mobs_f1_ + obs_f1_, 2)/std::pow(sigma_obs_f1_, 2);
            break;
        case 5:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_f2_, 2)) + 0.5*std::pow(-mobs_f2_ + obs_f2_, 2)/std::pow(sigma_obs_f2_, 2);
            break;
    }
}

} // namespace model_model_calvetti_py
} // namespace amici
