#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "z.h"
#include "sigmaz.h"
#include "mz.h"

namespace amici {
namespace model_model_neuron_py {

void dJzdsigma_model_neuron_py(realtype *dJzdsigma, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const realtype *mz){
    switch(iz) {
        case 0:
            dJzdsigma[0] = 1.0/sigma_z1 - 1.0*std::pow(-mz1 + z1, 2)/std::pow(sigma_z1, 3);
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
