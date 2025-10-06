#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"

namespace amici {
namespace model_model_neuron_py {

void dzdx_model_neuron_py(realtype *dzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h){
    switch(ie) {
        case 0:
            dzdx[0] = -1/(I0 - u + (1.0/25.0)*std::pow(v, 2) + 5*v + 140);
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
