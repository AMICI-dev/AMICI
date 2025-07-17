#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "k.h"

namespace amici {
namespace model_model_calvetti_py {

void x0_model_calvetti_py(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = V1ss;
    x0[1] = V2ss;
    x0[2] = V3ss;
    x0[3] = 1.0;
    x0[4] = 1.0;
    x0[5] = 1.0;
}

} // namespace model_model_calvetti_py
} // namespace amici
