#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_calvetti_py {

static constexpr std::array<sunindextype, 17> dwdw_colptrs_model_calvetti_py_ = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 14, 18, 18, 18, 18, 18
};

void dwdw_colptrs_model_calvetti_py(SUNMatrixWrapper &dwdw){
    dwdw.set_indexptrs(gsl::make_span(dwdw_colptrs_model_calvetti_py_));
}
} // namespace model_model_calvetti_py
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_calvetti_py {

static constexpr std::array<sunindextype, 18> dwdw_rowvals_model_calvetti_py_ = {
    13, 9, 10, 11, 7, 8, 13, 14, 15, 12, 13, 12, 13, 14, 12, 13, 14, 15
};

void dwdw_rowvals_model_calvetti_py(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_model_calvetti_py_));
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

void dwdw_model_calvetti_py(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static){
    const realtype V1_ = x[0];
    const realtype V2_ = x[1];
    const realtype V3_ = x[2];
    const realtype f1_ = x[3];
    const realtype f2_ = x[4];
    const realtype f3_ = x[5];
    const realtype V2ss_ = k[2];
    const realtype R2ss_ = k[3];
    const realtype V3ss_ = k[4];
    const realtype R3ss_ = k[5];
    const realtype C1ss_ = w[0];
    const realtype L1_ = w[1];
    const realtype L2_ = w[2];
    const realtype L3_ = w[3];
    const realtype p2_ = w[4];
    const realtype p3_ = w[5];
    const realtype C2ss_ = w[7];
    const realtype C3ss_ = w[8];
    const realtype R1_ = w[9];
    const realtype R2_ = w[10];
    const realtype R3_ = w[11];

    realtype &drate_of_V1_dC1ss_ = dwdw[0];
    realtype &dR1_dL1_ = dwdw[1];
    realtype &dR2_dL2_ = dwdw[2];
    realtype &dR3_dL3_ = dwdw[3];
    realtype &dC2ss_dp2_ = dwdw[4];
    realtype &dC3ss_dp3_ = dwdw[5];
    realtype &drate_of_V1_ds_ = dwdw[6];
    realtype &drate_of_V2_dC2ss_ = dwdw[7];
    realtype &drate_of_V3_dC3ss_ = dwdw[8];
    realtype &df0_dR1_ = dwdw[9];
    realtype &drate_of_V1_dR1_ = dwdw[10];
    realtype &df0_dR2_ = dwdw[11];
    realtype &drate_of_V1_dR2_ = dwdw[12];
    realtype &drate_of_V2_dR2_ = dwdw[13];
    realtype &df0_dR3_ = dwdw[14];
    realtype &drate_of_V1_dR3_ = dwdw[15];
    realtype &drate_of_V2_dR3_ = dwdw[16];
    realtype &drate_of_V3_dR3_ = dwdw[17];
    // static expressions
    if (include_static) {
        dC2ss_dp2_ = -V2ss_/std::pow(-0.5*R2ss_ + p2_, 2);  // dwdw[4]
        dC3ss_dp3_ = -V3ss_/std::pow(-0.5*R3ss_ + p3_, 2);  // dwdw[5]
        drate_of_V1_ds_ = 1.0/31.0;  // dwdw[6]
    }

    // dynamic expressions
    drate_of_V1_dC1ss_ = (2.0/31.0)*V1_/(std::pow(C1ss_, 2)*(R3_*f3_ + f1_*(R1_ + R2_) + f2_*(R2_ + R3_)));  // dwdw[0]
    dR1_dL1_ = 3*std::pow(L1_, 2)/std::pow(V1_, 2);  // dwdw[1]
    dR2_dL2_ = 3*std::pow(L2_, 2)/std::pow(V2_, 2);  // dwdw[2]
    dR3_dL3_ = 3*std::pow(L3_, 2)/std::pow(V3_, 2);  // dwdw[3]
    drate_of_V2_dC2ss_ = (2.0/163.0)*V2_/(std::pow(C2ss_, 2)*(R3_*f3_ + f2_*(R2_ + R3_)));  // dwdw[7]
    drate_of_V3_dC3ss_ = (1.0/61.0)*V3_/(std::pow(C3ss_, 2)*R3_*f3_);  // dwdw[8]
    df0_dR1_ = -f1_/R1_ + (R3_*f3_ + f1_*(R1_ + R2_) + f2_*(R2_ + R3_))/std::pow(R1_, 2) - 2/std::pow(R1_, 2);  // dwdw[9]
    drate_of_V1_dR1_ = (2.0/31.0)*V1_*f1_/(C1ss_*std::pow(R3_*f3_ + f1_*(R1_ + R2_) + f2_*(R2_ + R3_), 2));  // dwdw[10]
    df0_dR2_ = (-f1_ - f2_)/R1_;  // dwdw[11]
    drate_of_V1_dR2_ = (2.0/31.0)*V1_*(f1_ + f2_)/(C1ss_*std::pow(R3_*f3_ + f1_*(R1_ + R2_) + f2_*(R2_ + R3_), 2));  // dwdw[12]
    drate_of_V2_dR2_ = (2.0/163.0)*V2_*f2_/(C2ss_*std::pow(R3_*f3_ + f2_*(R2_ + R3_), 2));  // dwdw[13]
    df0_dR3_ = (-f2_ - f3_)/R1_;  // dwdw[14]
    drate_of_V1_dR3_ = (2.0/31.0)*V1_*(f2_ + f3_)/(C1ss_*std::pow(R3_*f3_ + f1_*(R1_ + R2_) + f2_*(R2_ + R3_), 2));  // dwdw[15]
    drate_of_V2_dR3_ = (2.0/163.0)*V2_*(f2_ + f3_)/(C2ss_*std::pow(R3_*f3_ + f2_*(R2_ + R3_), 2));  // dwdw[16]
    drate_of_V3_dR3_ = (1.0/61.0)*V3_/(C3ss_*std::pow(R3_, 2)*f3_);  // dwdw[17]
}

} // namespace model_model_calvetti_py
} // namespace amici
