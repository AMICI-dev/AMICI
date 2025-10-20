#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "p.h"
#include "k.h"
#include "y.h"
#include "sigmay.h"
#include "my.h"

namespace amici {
namespace model_model_steadystate_py {

void Jy_model_steadystate_py(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_x1, 2)) + 0.5*std::pow(-mobs_x1 + obs_x1, 2)/std::pow(sigma_obs_x1, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_x2, 2)) + 0.5*std::pow(-mobs_x2 + obs_x2, 2)/std::pow(sigma_obs_x2, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_x3, 2)) + 0.5*std::pow(-mobs_x3 + obs_x3, 2)/std::pow(sigma_obs_x3, 2);
            break;
    }
}

} // namespace model_model_steadystate_py
} // namespace amici
