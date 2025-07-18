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
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"
#include "spl.h"
#include "sspl.h"
#include "dwdp.h"

namespace amici {
namespace model_model_jakstat_adjoint_py {

void dwdp_model_jakstat_adjoint_py(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){

    // dynamic expressions
    du_dsp1 = sspl_0_5;  // dwdp[0]
    du_dsp2 = sspl_0_6;  // dwdp[1]
    du_dsp3 = sspl_0_7;  // dwdp[2]
    du_dsp4 = sspl_0_8;  // dwdp[3]
    du_dsp5 = sspl_0_9;  // dwdp[4]
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
