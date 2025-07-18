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
namespace model_model_events_py {

void sigmay_model_events_py(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_y1 = 1.0;  // sigmay[0]
}

} // namespace model_model_events_py
} // namespace amici
