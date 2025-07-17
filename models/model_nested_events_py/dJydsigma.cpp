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

void dJydsigma_model_nested_events_py(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_obs_Virus - 1.0*std::pow(-mobs_Virus + obs_Virus, 2)/std::pow(sigma_obs_Virus, 3);
            break;
    }
}

} // namespace model_model_nested_events_py
} // namespace amici
