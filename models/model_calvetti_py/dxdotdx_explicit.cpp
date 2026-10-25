#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_calvetti_py {

static constexpr std::array<sunindextype, 7> dxdotdx_explicit_colptrs_model_calvetti_py_ = {
    0, 0, 0, 0, 2, 4, 5
};

void dxdotdx_explicit_colptrs_model_calvetti_py(SUNMatrixWrapper &dxdotdx_explicit){
    dxdotdx_explicit.set_indexptrs(gsl::make_span(dxdotdx_explicit_colptrs_model_calvetti_py_));
}
} // namespace model_model_calvetti_py
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_calvetti_py {

static constexpr std::array<sunindextype, 5> dxdotdx_explicit_rowvals_model_calvetti_py_ = {
    3, 4, 4, 5, 5
};

void dxdotdx_explicit_rowvals_model_calvetti_py(SUNMatrixWrapper &dxdotdx_explicit){
    dxdotdx_explicit.set_indexvals(gsl::make_span(dxdotdx_explicit_rowvals_model_calvetti_py_));
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

void dxdotdx_explicit_model_calvetti_py(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *w){
    realtype &dae_0_df1_ = dxdotdx_explicit[0];
    realtype &dae_1_df1_ = dxdotdx_explicit[1];
    realtype &dae_1_df2_ = dxdotdx_explicit[2];
    realtype &dae_2_df2_ = dxdotdx_explicit[3];
    realtype &dae_2_df3_ = dxdotdx_explicit[4];
    dae_0_df1_ = -1;  // dxdotdx_explicit[0]
    dae_1_df1_ = 1;  // dxdotdx_explicit[1]
    dae_1_df2_ = -1;  // dxdotdx_explicit[2]
    dae_2_df2_ = 1;  // dxdotdx_explicit[3]
    dae_2_df3_ = -1;  // dxdotdx_explicit[4]
}

} // namespace model_model_calvetti_py
} // namespace amici
