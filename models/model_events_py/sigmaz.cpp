#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_events_py {

void sigmaz_model_events_py(realtype *sigmaz, const realtype t, const realtype *p, const realtype *k){
    realtype &sigma_z1_ = sigmaz[0];
    realtype &sigma_z2_ = sigmaz[1];
    sigma_z1_ = 1.0;  // sigmaz[0]
    sigma_z2_ = 1.0;  // sigmaz[1]
}

} // namespace model_model_events_py
} // namespace amici
