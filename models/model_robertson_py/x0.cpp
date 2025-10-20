#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "p.h"
#include "k.h"

namespace amici {
namespace model_model_robertson_py {

void x0_model_robertson_py(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = k1;
}

} // namespace model_model_robertson_py
} // namespace amici
