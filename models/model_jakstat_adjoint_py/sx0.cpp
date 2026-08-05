#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_model_jakstat_adjoint_py {

void sx0_model_jakstat_adjoint_py(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip){
    switch(ip) {
        case 4:
            sx0[0] = 1;
            break;
    }
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
