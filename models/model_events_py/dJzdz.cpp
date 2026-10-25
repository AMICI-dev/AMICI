#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_events_py {

void dJzdz_model_events_py(realtype *dJzdz, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const double *mz){
    const realtype z1_ = z[0];
    const realtype z2_ = z[1];
    const realtype sigma_z1_ = sigmaz[0];
    const realtype sigma_z2_ = sigmaz[1];
    const realtype mz1_ = mz[0];
    const realtype mz2_ = mz[1];

    switch(iz) {
        case 0:
            dJzdz[0] = (-1.0*mz1_ + 1.0*z1_)/std::pow(sigma_z1_, 2);
            break;
        case 1:
            dJzdz[1] = (-1.0*mz2_ + 1.0*z2_)/std::pow(sigma_z2_, 2);
            break;
    }
}

} // namespace model_model_events_py
} // namespace amici
