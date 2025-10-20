#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "p.h"
#include "y.h"
#include "sigmay.h"
#include "my.h"

namespace amici {
namespace model_model_dirac_py {

void dJydsigma_model_dirac_py(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_obs_x2 - 1.0*std::pow(-mobs_x2 + obs_x2, 2)/std::pow(sigma_obs_x2, 3);
            break;
    }
}

} // namespace model_model_dirac_py
} // namespace amici
