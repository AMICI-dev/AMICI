#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_calvetti_py {

static constexpr std::array<std::array<sunindextype, 7>, 6> dJydy_colptrs_model_calvetti_py_ = {{
    {0, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 1}, 
}};

void dJydy_colptrs_model_calvetti_py(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_model_calvetti_py_[index]));
}
} // namespace model_model_calvetti_py
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_calvetti_py {

static constexpr std::array<std::array<sunindextype, 1>, 6> dJydy_rowvals_model_calvetti_py_ = {{
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_model_calvetti_py(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_model_calvetti_py_[index]));
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

void dJydy_model_calvetti_py(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    const realtype obs_V1_ = y[0];
    const realtype obs_V2_ = y[1];
    const realtype obs_V3_ = y[2];
    const realtype obs_f0_ = y[3];
    const realtype obs_f1_ = y[4];
    const realtype obs_f2_ = y[5];
    const realtype sigma_obs_V1_ = sigmay[0];
    const realtype sigma_obs_V2_ = sigmay[1];
    const realtype sigma_obs_V3_ = sigmay[2];
    const realtype sigma_obs_f0_ = sigmay[3];
    const realtype sigma_obs_f1_ = sigmay[4];
    const realtype sigma_obs_f2_ = sigmay[5];
    const realtype mobs_V1_ = my[0];
    const realtype mobs_V2_ = my[1];
    const realtype mobs_V3_ = my[2];
    const realtype mobs_f0_ = my[3];
    const realtype mobs_f1_ = my[4];
    const realtype mobs_f2_ = my[5];

    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mobs_V1_ + 1.0*obs_V1_)/std::pow(sigma_obs_V1_, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*mobs_V2_ + 1.0*obs_V2_)/std::pow(sigma_obs_V2_, 2);
            break;
        case 2:
            dJydy[0] = (-1.0*mobs_V3_ + 1.0*obs_V3_)/std::pow(sigma_obs_V3_, 2);
            break;
        case 3:
            dJydy[0] = (-1.0*mobs_f0_ + 1.0*obs_f0_)/std::pow(sigma_obs_f0_, 2);
            break;
        case 4:
            dJydy[0] = (-1.0*mobs_f1_ + 1.0*obs_f1_)/std::pow(sigma_obs_f1_, 2);
            break;
        case 5:
            dJydy[0] = (-1.0*mobs_f2_ + 1.0*obs_f2_)/std::pow(sigma_obs_f2_, 2);
            break;
    }
}

} // namespace model_model_calvetti_py
} // namespace amici
