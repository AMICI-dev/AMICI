#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_jakstat_adjoint_py {

static constexpr std::array<sunindextype, 18> dwdp_colptrs_model_jakstat_adjoint_py_ = {
    0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5, 5, 5
};

void dwdp_colptrs_model_jakstat_adjoint_py(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_model_jakstat_adjoint_py_));
}
} // namespace model_model_jakstat_adjoint_py
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_jakstat_adjoint_py {

static constexpr std::array<sunindextype, 5> dwdp_rowvals_model_jakstat_adjoint_py_ = {
    0, 0, 0, 0, 0
};

void dwdp_rowvals_model_jakstat_adjoint_py(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_model_jakstat_adjoint_py_));
}
} // namespace model_model_jakstat_adjoint_py
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include <sundials/sundials_types.h>
#include <gsl/gsl-lite.hpp>

namespace amici {
namespace model_model_jakstat_adjoint_py {

void dwdp_model_jakstat_adjoint_py(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){
    const realtype sspl_0_5_ = sspl[5];
    const realtype sspl_0_6_ = sspl[6];
    const realtype sspl_0_7_ = sspl[7];
    const realtype sspl_0_8_ = sspl[8];
    const realtype sspl_0_9_ = sspl[9];

    realtype &du_dsp1_ = dwdp[0];
    realtype &du_dsp2_ = dwdp[1];
    realtype &du_dsp3_ = dwdp[2];
    realtype &du_dsp4_ = dwdp[3];
    realtype &du_dsp5_ = dwdp[4];

    // dynamic expressions
    du_dsp1_ = sspl_0_5_;  // dwdp[0]
    du_dsp2_ = sspl_0_6_;  // dwdp[1]
    du_dsp3_ = sspl_0_7_;  // dwdp[2]
    du_dsp4_ = sspl_0_8_;  // dwdp[3]
    du_dsp5_ = sspl_0_9_;  // dwdp[4]
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
