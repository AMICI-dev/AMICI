#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_dirac_py {

void dJydsigmay_model_dirac_py(realtype *dJydsigmay, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    const realtype obs_x2_ = y[0];
    const realtype sigma_obs_x2_ = sigmay[0];
    const realtype mobs_x2_ = my[0];

    switch(iy) {
        case 0:
            dJydsigmay[0] = 1.0/sigma_obs_x2_ - 1.0*std::pow(-mobs_x2_ + obs_x2_, 2)/std::pow(sigma_obs_x2_, 3);
            break;
    }
}

} // namespace model_model_dirac_py
} // namespace amici
