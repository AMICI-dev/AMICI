#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "p.h"
#include "k.h"
#include "z.h"
#include "sigmaz.h"
#include "mz.h"

namespace amici {
namespace model_model_neuron_py {

void dJzdz_model_neuron_py(realtype *dJzdz, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const double *mz){
    switch(iz) {
        case 0:
            dJzdz[0] = (-1.0*mz1 + 1.0*z1)/std::pow(sigma_z1, 2);
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
