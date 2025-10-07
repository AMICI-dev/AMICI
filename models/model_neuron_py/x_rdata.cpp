#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_model_neuron_py {

void x_rdata_model_neuron_py(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = v;
    x_rdata[1] = u;
}

} // namespace model_model_neuron_py
} // namespace amici
