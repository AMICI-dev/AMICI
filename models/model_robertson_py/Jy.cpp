#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_robertson_py {

void Jy_model_robertson_py(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    const realtype obs_x1_ = y[0];
    const realtype obs_x2_ = y[1];
    const realtype obs_x3_ = y[2];
    const realtype sigma_obs_x1_ = sigmay[0];
    const realtype sigma_obs_x2_ = sigmay[1];
    const realtype sigma_obs_x3_ = sigmay[2];
    const realtype mobs_x1_ = my[0];
    const realtype mobs_x2_ = my[1];
    const realtype mobs_x3_ = my[2];

    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_x1_, 2)) + 0.5*std::pow(-mobs_x1_ + obs_x1_, 2)/std::pow(sigma_obs_x1_, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_x2_, 2)) + 0.5*std::pow(-mobs_x2_ + obs_x2_, 2)/std::pow(sigma_obs_x2_, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_x3_, 2)) + 0.5*std::pow(-mobs_x3_ + obs_x3_, 2)/std::pow(sigma_obs_x3_, 2);
            break;
    }
}

} // namespace model_model_robertson_py
} // namespace amici
