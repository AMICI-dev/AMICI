#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

void dzdt_model_neuron_py(realtype *dzdt, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h){
    switch(ie) {
        case 0:
            dzdt[0] = 1;
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
