#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"
#include "z.h"

namespace amici {
namespace model_model_neuron_py {

void z_model_neuron_py(realtype *z, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h){
    switch(ie) {
        case 0:
            z[0] = t;
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
