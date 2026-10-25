#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_jakstat_adjoint_py {

void dydp_model_jakstat_adjoint_py(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl){
    const realtype STAT_ = x[0];
    const realtype pSTAT_ = x[1];
    const realtype pSTAT_pSTAT_ = x[2];
    const realtype init_STAT_ = p[4];
    const realtype scale_tSTAT_ = p[12];
    const realtype scale_pSTAT_ = p[13];
    const realtype sspl_0_5_ = sspl[5];
    const realtype sspl_0_6_ = sspl[6];
    const realtype sspl_0_7_ = sspl[7];
    const realtype sspl_0_8_ = sspl[8];
    const realtype sspl_0_9_ = sspl[9];

    switch(ip) {
        case 4:
            dydp[0] = -scale_pSTAT_*(pSTAT_ + 2*pSTAT_pSTAT_)/std::pow(init_STAT_, 2);
            dydp[1] = -scale_tSTAT_*(STAT_ + pSTAT_ + 2*pSTAT_pSTAT_)/std::pow(init_STAT_, 2);
            break;
        case 5:
            dydp[2] = sspl_0_5_;
            break;
        case 6:
            dydp[2] = sspl_0_6_;
            break;
        case 7:
            dydp[2] = sspl_0_7_;
            break;
        case 8:
            dydp[2] = sspl_0_8_;
            break;
        case 9:
            dydp[2] = sspl_0_9_;
            break;
        case 10:
            dydp[1] = 1;
            break;
        case 11:
            dydp[0] = 1;
            break;
        case 12:
            dydp[1] = (STAT_ + pSTAT_ + 2*pSTAT_pSTAT_)/init_STAT_;
            break;
        case 13:
            dydp[0] = (pSTAT_ + 2*pSTAT_pSTAT_)/init_STAT_;
            break;
    }
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
