#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_jakstat_adjoint_py {

void y_model_jakstat_adjoint_py(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype STAT_ = x[0];
    const realtype pSTAT_ = x[1];
    const realtype pSTAT_pSTAT_ = x[2];
    const realtype init_STAT_ = p[4];
    const realtype offset_tSTAT_ = p[10];
    const realtype offset_pSTAT_ = p[11];
    const realtype scale_tSTAT_ = p[12];
    const realtype scale_pSTAT_ = p[13];
    const realtype u_ = w[0];

    y[0] = offset_pSTAT_ + scale_pSTAT_*(pSTAT_ + 2*pSTAT_pSTAT_)/init_STAT_;
    y[1] = offset_tSTAT_ + scale_tSTAT_*(STAT_ + pSTAT_ + 2*pSTAT_pSTAT_)/init_STAT_;
    y[2] = u_;
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
