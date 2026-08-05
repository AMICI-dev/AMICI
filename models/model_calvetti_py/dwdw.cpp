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
#include "x.h"
#include "k.h"
#include "h.h"
#include "w.h"
#include "dwdw.h"

namespace amici {
namespace model_model_calvetti_py {

void dwdw_model_calvetti_py(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static){
    // static expressions
    if (include_static) {
        dC2ss_dp2 = -V2ss/std::pow(-0.5*R2ss + p2, 2);  // dwdw[4]
        dC3ss_dp3 = -V3ss/std::pow(-0.5*R3ss + p3, 2);  // dwdw[5]
        drate_of_V1_ds = 1.0/31.0;  // dwdw[6]
    }

    // dynamic expressions
    drate_of_V1_dC1ss = (2.0/31.0)*V1/(std::pow(C1ss, 2)*(R3*f3 + f1*(R1 + R2) + f2*(R2 + R3)));  // dwdw[0]
    dR1_dL1 = 3*std::pow(L1, 2)/std::pow(V1, 2);  // dwdw[1]
    dR2_dL2 = 3*std::pow(L2, 2)/std::pow(V2, 2);  // dwdw[2]
    dR3_dL3 = 3*std::pow(L3, 2)/std::pow(V3, 2);  // dwdw[3]
    drate_of_V2_dC2ss = (2.0/163.0)*V2/(std::pow(C2ss, 2)*(R3*f3 + f2*(R2 + R3)));  // dwdw[7]
    drate_of_V3_dC3ss = (1.0/61.0)*V3/(std::pow(C3ss, 2)*R3*f3);  // dwdw[8]
    df0_dR1 = -f1/R1 + (R3*f3 + f1*(R1 + R2) + f2*(R2 + R3))/std::pow(R1, 2) - 2/std::pow(R1, 2);  // dwdw[9]
    drate_of_V1_dR1 = (2.0/31.0)*V1*f1/(C1ss*std::pow(R3*f3 + f1*(R1 + R2) + f2*(R2 + R3), 2));  // dwdw[10]
    df0_dR2 = (-f1 - f2)/R1;  // dwdw[11]
    drate_of_V1_dR2 = (2.0/31.0)*V1*(f1 + f2)/(C1ss*std::pow(R3*f3 + f1*(R1 + R2) + f2*(R2 + R3), 2));  // dwdw[12]
    drate_of_V2_dR2 = (2.0/163.0)*V2*f2/(C2ss*std::pow(R3*f3 + f2*(R2 + R3), 2));  // dwdw[13]
    df0_dR3 = (-f2 - f3)/R1;  // dwdw[14]
    drate_of_V1_dR3 = (2.0/31.0)*V1*(f2 + f3)/(C1ss*std::pow(R3*f3 + f1*(R1 + R2) + f2*(R2 + R3), 2));  // dwdw[15]
    drate_of_V2_dR3 = (2.0/163.0)*V2*(f2 + f3)/(C2ss*std::pow(R3*f3 + f2*(R2 + R3), 2));  // dwdw[16]
    drate_of_V3_dR3 = (1.0/61.0)*V3/(C3ss*std::pow(R3, 2)*f3);  // dwdw[17]
}

} // namespace model_model_calvetti_py
} // namespace amici
