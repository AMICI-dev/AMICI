#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "p.h"
#include "k.h"
#include "sigmaz.h"

namespace amici {
namespace model_model_events_py {

void sigmaz_model_events_py(realtype *sigmaz, const realtype t, const realtype *p, const realtype *k){
    sigma_z1 = 1.0;  // sigmaz[0]
    sigma_z2 = 1.0;  // sigmaz[1]
}

} // namespace model_model_events_py
} // namespace amici
