#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_nested_events_py {

void x0_model_nested_events_py(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    const realtype V_0_ = p[0];

    x0[0] = V_0_;
}

} // namespace model_model_nested_events_py
} // namespace amici
