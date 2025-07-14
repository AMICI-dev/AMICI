#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_model_jakstat_adjoint_py {

void dspline_valuesdp_model_jakstat_adjoint_py(realtype *dspline_valuesdp, const realtype *p, const realtype *k, const int ip){
    switch(ip) {
        case 5:
            dspline_valuesdp[0] = 1;
            break;
        case 6:
            dspline_valuesdp[1] = 1;
            break;
        case 7:
            dspline_valuesdp[2] = 1;
            break;
        case 8:
            dspline_valuesdp[3] = 1;
            break;
        case 9:
            dspline_valuesdp[4] = 1;
            break;
    }
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
