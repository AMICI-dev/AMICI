#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"
#include "sigmay.h"

namespace amici {
namespace model_model_steadystate_py {

void sigmay_model_steadystate_py(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_obs_x1 = 1.0;  // sigmay[0]
    sigma_obs_x2 = 1.0;  // sigmay[1]
    sigma_obs_x3 = 1.0;  // sigmay[2]
}

} // namespace model_model_steadystate_py
} // namespace amici
