#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "y.h"
#include "sigmay.h"
#include "my.h"

namespace amici {
namespace model_model_nested_events_py {

void Jy_model_nested_events_py(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obs_Virus, 2)) + 0.5*std::pow(-mobs_Virus + obs_Virus, 2)/std::pow(sigma_obs_Virus, 2);
            break;
    }
}

} // namespace model_model_nested_events_py
} // namespace amici
