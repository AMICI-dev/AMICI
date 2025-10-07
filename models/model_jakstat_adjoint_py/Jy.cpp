#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "p.h"
#include "k.h"
#include "y.h"
#include "sigmay.h"
#include "my.h"

namespace amici {
namespace model_model_jakstat_adjoint_py {

void Jy_model_jakstat_adjoint_py(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_pSTAT, 2)) + 0.5*std::pow(-mobs_pSTAT + obs_pSTAT, 2)/std::pow(sigma_obs_pSTAT, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_tSTAT, 2)) + 0.5*std::pow(-mobs_tSTAT + obs_tSTAT, 2)/std::pow(sigma_obs_tSTAT, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_spline, 2)) + 0.5*std::pow(-mobs_spline + obs_spline, 2)/std::pow(sigma_obs_spline, 2);
            break;
    }
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
