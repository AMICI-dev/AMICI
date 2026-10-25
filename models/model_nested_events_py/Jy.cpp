#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_nested_events_py {

void Jy_model_nested_events_py(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    const realtype obs_Virus_ = y[0];
    const realtype sigma_obs_Virus_ = sigmay[0];
    const realtype mobs_Virus_ = my[0];

    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_Virus_, 2)) + 0.5*std::pow(-mobs_Virus_ + obs_Virus_, 2)/std::pow(sigma_obs_Virus_, 2);
            break;
    }
}

} // namespace model_model_nested_events_py
} // namespace amici
