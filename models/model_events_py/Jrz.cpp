#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "rz.h"
#include "sigmaz.h"

namespace amici {
namespace model_model_events_py {

void Jrz_model_events_py(realtype *Jrz, const int iz, const realtype *p, const realtype *k, const realtype *rz, const realtype *sigmaz){
    switch(iz) {
        case 0:
            Jrz[0] = 0.5*std::pow(rz1, 2)/std::pow(sigma_z1, 2) + 0.5*std::log(2*amici::pi*std::pow(sigma_z1, 2));
            Jrz[1] = 0.5*std::pow(rz2, 2)/std::pow(sigma_z2, 2) + 0.5*std::log(2*amici::pi*std::pow(sigma_z2, 2));
            break;
    }
}

} // namespace model_model_events_py
} // namespace amici
