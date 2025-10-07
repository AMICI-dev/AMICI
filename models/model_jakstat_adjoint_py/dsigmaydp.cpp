#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "p.h"
#include "k.h"
#include "y.h"

namespace amici {
namespace model_model_jakstat_adjoint_py {

void dsigmaydp_model_jakstat_adjoint_py(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 14:
            dsigmaydp[0] = 1;
            break;
        case 15:
            dsigmaydp[1] = 1;
            break;
        case 16:
            dsigmaydp[2] = 1;
            break;
    }
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
