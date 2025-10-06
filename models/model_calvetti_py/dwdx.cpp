#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_calvetti_py {

static constexpr std::array<sunindextype, 7> dwdx_colptrs_model_calvetti_py_ = {
    0, 2, 4, 6, 8, 11, 15
};

void dwdx_colptrs_model_calvetti_py(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_model_calvetti_py_));
}
} // namespace model_model_calvetti_py
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_calvetti_py {

static constexpr std::array<sunindextype, 15> dwdx_rowvals_model_calvetti_py_ = {
    9, 13, 10, 14, 11, 15, 12, 13, 12, 13, 14, 12, 13, 14, 15
};

void dwdx_rowvals_model_calvetti_py(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_model_calvetti_py_));
}
} // namespace model_model_calvetti_py
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include <sundials/sundials_types.h>
#include <gsl/gsl-lite.hpp>
#include "x.h"
#include "k.h"
#include "h.h"
#include "w.h"
#include "dwdx.h"

namespace amici {
namespace model_model_calvetti_py {

void dwdx_model_calvetti_py(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){

    // dynamic expressions
    dR1_dV1 = -2*std::pow(L1, 3)/std::pow(V1, 3);  // dwdx[0]
    drate_of_V1_dV1 = -(100.0/899.0)/V1ss - (2.0/31.0)/(C1ss*(R3*f3 + f1*(R1 + R2) + f2*(R2 + R3)));  // dwdx[1]
    dR2_dV2 = -2*std::pow(L2, 3)/std::pow(V2, 3);  // dwdx[2]
    drate_of_V2_dV2 = -(100.0/8313.0)/V2ss - (2.0/163.0)/(C2ss*(R3*f3 + f2*(R2 + R3)));  // dwdx[3]
    dR3_dV3 = -2*std::pow(L3, 3)/std::pow(V3, 3);  // dwdx[4]
    drate_of_V3_dV3 = -(1.0/121999878.0)/V3ss - (1.0/61.0)/(C3ss*R3*f3);  // dwdx[5]
    df0_df1 = (-R1 - R2)/R1;  // dwdx[6]
    drate_of_V1_df1 = (2.0/31.0)*V1*(R1 + R2)/(C1ss*std::pow(R3*f3 + f1*(R1 + R2) + f2*(R2 + R3), 2));  // dwdx[7]
    df0_df2 = (-R2 - R3)/R1;  // dwdx[8]
    drate_of_V1_df2 = (2.0/31.0)*V1*(R2 + R3)/(C1ss*std::pow(R3*f3 + f1*(R1 + R2) + f2*(R2 + R3), 2));  // dwdx[9]
    drate_of_V2_df2 = (2.0/163.0)*V2*(R2 + R3)/(C2ss*std::pow(R3*f3 + f2*(R2 + R3), 2));  // dwdx[10]
    df0_df3 = -R3/R1;  // dwdx[11]
    drate_of_V1_df3 = (2.0/31.0)*R3*V1/(C1ss*std::pow(R3*f3 + f1*(R1 + R2) + f2*(R2 + R3), 2));  // dwdx[12]
    drate_of_V2_df3 = (2.0/163.0)*R3*V2/(C2ss*std::pow(R3*f3 + f2*(R2 + R3), 2));  // dwdx[13]
    drate_of_V3_df3 = (1.0/61.0)*V3/(C3ss*R3*std::pow(f3, 2));  // dwdx[14]
}

} // namespace model_model_calvetti_py
} // namespace amici
