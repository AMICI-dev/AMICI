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
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "k.h"
#include "y.h"
#include "sigmay.h"
#include "my.h"
#include "dJydy.h"

namespace amici {
namespace model_model_calvetti_py {

void dJydy_model_calvetti_py(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mobs_V1 + 1.0*obs_V1)/std::pow(sigma_obs_V1, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*mobs_V2 + 1.0*obs_V2)/std::pow(sigma_obs_V2, 2);
            break;
        case 2:
            dJydy[0] = (-1.0*mobs_V3 + 1.0*obs_V3)/std::pow(sigma_obs_V3, 2);
            break;
        case 3:
            dJydy[0] = (-1.0*mobs_f0 + 1.0*obs_f0)/std::pow(sigma_obs_f0, 2);
            break;
        case 4:
            dJydy[0] = (-1.0*mobs_f1 + 1.0*obs_f1)/std::pow(sigma_obs_f1, 2);
            break;
        case 5:
            dJydy[0] = (-1.0*mobs_f2 + 1.0*obs_f2)/std::pow(sigma_obs_f2, 2);
            break;
    }
}

} // namespace model_model_calvetti_py
} // namespace amici
