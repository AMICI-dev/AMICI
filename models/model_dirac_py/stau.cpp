#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "h.h"
#include "sx.h"

namespace amici {
namespace model_model_dirac_py {

void stau_model_dirac_py(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *tcl, const realtype *sx, const int ip, const int ie){
    switch(ie) {
        case 0:
            switch(ip) {
                case 1:
                    stau[0] = -1;
                    break;
            }
            break;
    }
}

} // namespace model_model_dirac_py
} // namespace amici
