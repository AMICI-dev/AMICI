#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"
#include "sigmay.h"
#include "my.h"

namespace amici {
namespace model_model_robertson_py {

void dJydsigma_model_robertson_py(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_obs_x1 - 1.0*std::pow(-mobs_x1 + obs_x1, 2)/std::pow(sigma_obs_x1, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_obs_x2 - 1.0*std::pow(-mobs_x2 + obs_x2, 2)/std::pow(sigma_obs_x2, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_obs_x3 - 1.0*std::pow(-mobs_x3 + obs_x3, 2)/std::pow(sigma_obs_x3, 3);
            break;
    }
}

} // namespace model_model_robertson_py
} // namespace amici
