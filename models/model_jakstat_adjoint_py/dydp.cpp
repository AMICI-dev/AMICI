#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"
#include "spl.h"
#include "sspl.h"

namespace amici {
namespace model_model_jakstat_adjoint_py {

void dydp_model_jakstat_adjoint_py(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl){
    switch(ip) {
        case 4:
            dydp[0] = -scale_pSTAT*(pSTAT + 2*pSTAT_pSTAT)/std::pow(init_STAT, 2);
            dydp[1] = -scale_tSTAT*(STAT + pSTAT + 2*pSTAT_pSTAT)/std::pow(init_STAT, 2);
            break;
        case 5:
            dydp[2] = sspl_0_5;
            break;
        case 6:
            dydp[2] = sspl_0_6;
            break;
        case 7:
            dydp[2] = sspl_0_7;
            break;
        case 8:
            dydp[2] = sspl_0_8;
            break;
        case 9:
            dydp[2] = sspl_0_9;
            break;
        case 10:
            dydp[1] = 1;
            break;
        case 11:
            dydp[0] = 1;
            break;
        case 12:
            dydp[1] = (STAT + pSTAT + 2*pSTAT_pSTAT)/init_STAT;
            break;
        case 13:
            dydp[0] = (pSTAT + 2*pSTAT_pSTAT)/init_STAT;
            break;
    }
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
