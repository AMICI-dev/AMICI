#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_dirac_py {

void Jy_model_dirac_py(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    const realtype obs_x2_ = y[0];
    const realtype sigma_obs_x2_ = sigmay[0];
    const realtype mobs_x2_ = my[0];

    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_x2_, 2)) + 0.5*std::pow(-mobs_x2_ + obs_x2_, 2)/std::pow(sigma_obs_x2_, 2);
            break;
    }
}

} // namespace model_model_dirac_py
} // namespace amici
