#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "p.h"
#include "k.h"
#include "rz.h"
#include "sigmaz.h"

namespace amici {
namespace model_model_events_py {

void dJrzdz_model_events_py(realtype *dJrzdz, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz){
    switch(iz) {
        case 0:
            dJrzdz[0] = 1.0*rz1/std::pow(sigma_z1, 2);
            break;
        case 1:
            dJrzdz[1] = 1.0*rz2/std::pow(sigma_z2, 2);
            break;
    }
}

} // namespace model_model_events_py
} // namespace amici
