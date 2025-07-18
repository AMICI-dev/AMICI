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
namespace model_model_jakstat_adjoint_py {

void dJydsigma_model_jakstat_adjoint_py(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_obs_pSTAT - 1.0*std::pow(-mobs_pSTAT + obs_pSTAT, 2)/std::pow(sigma_obs_pSTAT, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_obs_tSTAT - 1.0*std::pow(-mobs_tSTAT + obs_tSTAT, 2)/std::pow(sigma_obs_tSTAT, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_obs_spline - 1.0*std::pow(-mobs_spline + obs_spline, 2)/std::pow(sigma_obs_spline, 3);
            break;
    }
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
