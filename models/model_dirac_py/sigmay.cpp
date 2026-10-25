#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_dirac_py {

void sigmay_model_dirac_py(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    realtype &sigma_obs_x2_ = sigmay[0];
    sigma_obs_x2_ = 1.0;  // sigmay[0]
}

} // namespace model_model_dirac_py
} // namespace amici
