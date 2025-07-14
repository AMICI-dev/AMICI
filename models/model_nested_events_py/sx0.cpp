#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"

namespace amici {
namespace model_model_nested_events_py {

void sx0_model_nested_events_py(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip){
    switch(ip) {
        case 0:
            sx0[0] = 1;
            break;
    }
}

} // namespace model_model_nested_events_py
} // namespace amici
