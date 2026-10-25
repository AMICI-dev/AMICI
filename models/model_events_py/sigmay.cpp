#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_events_py {

void sigmay_model_events_py(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    realtype &sigma_y1_ = sigmay[0];
    sigma_y1_ = 1.0;  // sigmay[0]
}

} // namespace model_model_events_py
} // namespace amici
