#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"

namespace amici {
namespace model_model_neuron_py {

void drzdx_model_neuron_py(realtype *drzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h){
    switch(ie) {
        case 0:
            drzdx[0] = 1;
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
