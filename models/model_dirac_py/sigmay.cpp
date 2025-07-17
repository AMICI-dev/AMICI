#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "y.h"
#include "sigmay.h"

namespace amici {
namespace model_model_dirac_py {

void sigmay_model_dirac_py(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_obs_x2 = 1.0;  // sigmay[0]
}

} // namespace model_model_dirac_py
} // namespace amici
