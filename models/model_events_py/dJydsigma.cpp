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
namespace model_model_events_py {

void dJydsigma_model_events_py(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_y1 - 1.0*std::pow(-my1 + y1, 2)/std::pow(sigma_y1, 3);
            break;
    }
}

} // namespace model_model_events_py
} // namespace amici
