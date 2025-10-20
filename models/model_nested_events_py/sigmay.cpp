#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "p.h"
#include "y.h"
#include "sigmay.h"

namespace amici {
namespace model_model_nested_events_py {

void sigmay_model_nested_events_py(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_obs_Virus = 1.0;  // sigmay[0]
}

} // namespace model_model_nested_events_py
} // namespace amici
