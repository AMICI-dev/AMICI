#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

void drootdx_model_neuron_py(realtype *drootdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const realtype *tcl, const realtype *spl){
    drootdx[0] = 1;
}

} // namespace model_model_neuron_py
} // namespace amici
