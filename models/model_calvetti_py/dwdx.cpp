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

namespace amici {
namespace model_model_calvetti_py {

void dwdx_model_calvetti_py(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    const realtype V1_ = x[0];
    const realtype V2_ = x[1];
    const realtype V3_ = x[2];
    const realtype f1_ = x[3];
    const realtype f2_ = x[4];
    const realtype f3_ = x[5];
    const realtype V1ss_ = k[0];
    const realtype V2ss_ = k[2];
    const realtype V3ss_ = k[4];
    const realtype C1ss_ = w[0];
    const realtype L1_ = w[1];
    const realtype L2_ = w[2];
    const realtype L3_ = w[3];
    const realtype C2ss_ = w[7];
    const realtype C3ss_ = w[8];
    const realtype R1_ = w[9];
    const realtype R2_ = w[10];
    const realtype R3_ = w[11];

    realtype &dR1_dV1_ = dwdx[0];
    realtype &drate_of_V1_dV1_ = dwdx[1];
    realtype &dR2_dV2_ = dwdx[2];
    realtype &drate_of_V2_dV2_ = dwdx[3];
    realtype &dR3_dV3_ = dwdx[4];
    realtype &drate_of_V3_dV3_ = dwdx[5];
    realtype &df0_df1_ = dwdx[6];
    realtype &drate_of_V1_df1_ = dwdx[7];
    realtype &df0_df2_ = dwdx[8];
    realtype &drate_of_V1_df2_ = dwdx[9];
    realtype &drate_of_V2_df2_ = dwdx[10];
    realtype &df0_df3_ = dwdx[11];
    realtype &drate_of_V1_df3_ = dwdx[12];
    realtype &drate_of_V2_df3_ = dwdx[13];
    realtype &drate_of_V3_df3_ = dwdx[14];

    // dynamic expressions
    dR1_dV1_ = -2*std::pow(L1_, 3)/std::pow(V1_, 3);  // dwdx[0]
    drate_of_V1_dV1_ = -(100.0/899.0)/V1ss_ - (2.0/31.0)/(C1ss_*(R3_*f3_ + f1_*(R1_ + R2_) + f2_*(R2_ + R3_)));  // dwdx[1]
    dR2_dV2_ = -2*std::pow(L2_, 3)/std::pow(V2_, 3);  // dwdx[2]
    drate_of_V2_dV2_ = -(100.0/8313.0)/V2ss_ - (2.0/163.0)/(C2ss_*(R3_*f3_ + f2_*(R2_ + R3_)));  // dwdx[3]
    dR3_dV3_ = -2*std::pow(L3_, 3)/std::pow(V3_, 3);  // dwdx[4]
    drate_of_V3_dV3_ = -(1.0/121999878.0)/V3ss_ - (1.0/61.0)/(C3ss_*R3_*f3_);  // dwdx[5]
    df0_df1_ = (-R1_ - R2_)/R1_;  // dwdx[6]
    drate_of_V1_df1_ = (2.0/31.0)*V1_*(R1_ + R2_)/(C1ss_*std::pow(R3_*f3_ + f1_*(R1_ + R2_) + f2_*(R2_ + R3_), 2));  // dwdx[7]
    df0_df2_ = (-R2_ - R3_)/R1_;  // dwdx[8]
    drate_of_V1_df2_ = (2.0/31.0)*V1_*(R2_ + R3_)/(C1ss_*std::pow(R3_*f3_ + f1_*(R1_ + R2_) + f2_*(R2_ + R3_), 2));  // dwdx[9]
    drate_of_V2_df2_ = (2.0/163.0)*V2_*(R2_ + R3_)/(C2ss_*std::pow(R3_*f3_ + f2_*(R2_ + R3_), 2));  // dwdx[10]
    df0_df3_ = -R3_/R1_;  // dwdx[11]
    drate_of_V1_df3_ = (2.0/31.0)*R3_*V1_/(C1ss_*std::pow(R3_*f3_ + f1_*(R1_ + R2_) + f2_*(R2_ + R3_), 2));  // dwdx[12]
    drate_of_V2_df3_ = (2.0/163.0)*R3_*V2_/(C2ss_*std::pow(R3_*f3_ + f2_*(R2_ + R3_), 2));  // dwdx[13]
    drate_of_V3_df3_ = (1.0/61.0)*V3_/(C3ss_*R3_*std::pow(f3_, 2));  // dwdx[14]
}

} // namespace model_model_calvetti_py
} // namespace amici
