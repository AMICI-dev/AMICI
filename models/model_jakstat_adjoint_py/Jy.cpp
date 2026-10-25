#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_jakstat_adjoint_py {

void Jy_model_jakstat_adjoint_py(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    const realtype obs_pSTAT_ = y[0];
    const realtype obs_tSTAT_ = y[1];
    const realtype obs_spline_ = y[2];
    const realtype sigma_obs_pSTAT_ = sigmay[0];
    const realtype sigma_obs_tSTAT_ = sigmay[1];
    const realtype sigma_obs_spline_ = sigmay[2];
    const realtype mobs_pSTAT_ = my[0];
    const realtype mobs_tSTAT_ = my[1];
    const realtype mobs_spline_ = my[2];

    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_pSTAT_, 2)) + 0.5*std::pow(-mobs_pSTAT_ + obs_pSTAT_, 2)/std::pow(sigma_obs_pSTAT_, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_tSTAT_, 2)) + 0.5*std::pow(-mobs_tSTAT_ + obs_tSTAT_, 2)/std::pow(sigma_obs_tSTAT_, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_spline_, 2)) + 0.5*std::pow(-mobs_spline_ + obs_spline_, 2)/std::pow(sigma_obs_spline_, 2);
            break;
    }
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
