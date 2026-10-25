#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_calvetti_py {

static constexpr std::array<sunindextype, 17> dxdotdw_colptrs_model_calvetti_py_ = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 5, 7
};

void dxdotdw_colptrs_model_calvetti_py(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_model_calvetti_py_));
}
} // namespace model_model_calvetti_py
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_calvetti_py {

static constexpr std::array<sunindextype, 7> dxdotdw_rowvals_model_calvetti_py_ = {
    3, 0, 3, 1, 4, 2, 5
};

void dxdotdw_rowvals_model_calvetti_py(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_model_calvetti_py_));
}
} // namespace model_model_calvetti_py
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include <sundials/sundials_types.h>
#include <gsl/gsl-lite.hpp>

namespace amici {
namespace model_model_calvetti_py {

void dxdotdw_model_calvetti_py(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *w){
    realtype &dae_0_df0_ = dxdotdw[0];
    realtype &dde_0_drate_of_V1_ = dxdotdw[1];
    realtype &dae_0_drate_of_V1_ = dxdotdw[2];
    realtype &dde_1_drate_of_V2_ = dxdotdw[3];
    realtype &dae_1_drate_of_V2_ = dxdotdw[4];
    realtype &dde_2_drate_of_V3_ = dxdotdw[5];
    realtype &dae_2_drate_of_V3_ = dxdotdw[6];
    dae_0_df0_ = 1;  // dxdotdw[0]
    dde_0_drate_of_V1_ = 1;  // dxdotdw[1]
    dae_0_drate_of_V1_ = -1;  // dxdotdw[2]
    dde_1_drate_of_V2_ = 1;  // dxdotdw[3]
    dae_1_drate_of_V2_ = -1;  // dxdotdw[4]
    dde_2_drate_of_V3_ = 1;  // dxdotdw[5]
    dae_2_drate_of_V3_ = -1;  // dxdotdw[6]
}

} // namespace model_model_calvetti_py
} // namespace amici
