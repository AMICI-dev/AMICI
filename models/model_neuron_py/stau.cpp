#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"
#include "sx.h"

namespace amici {
namespace model_model_neuron_py {

void stau_model_neuron_py(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const realtype *tcl, const realtype *sx, const int ip, const int ie){
    switch(ie) {
        case 0:
            switch(ip) {
                case 0:
                case 1:
                case 2:
                case 3:
                    stau[0] = sx0/(I0 - u + (1.0/25.0)*std::pow(v, 2) + 5*v + 140);
                    break;
            }
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
