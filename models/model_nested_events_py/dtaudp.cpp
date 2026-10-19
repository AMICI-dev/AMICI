#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_nested_events_py {

void dtaudp_model_nested_events_py(realtype *dtaudp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dx, const realtype *tcl, const int ip, const int ie){
    switch(ie) {
        case 2:
            switch(ip) {
                case 2:
                    dtaudp[0] = -1;
                    break;
            }
            break;
    }
}

} // namespace model_model_nested_events_py
} // namespace amici
