#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_robertson_py {

static constexpr std::array<std::array<sunindextype, 4>, 3> dJydy_colptrs_model_robertson_py_ = {{
    {0, 1, 1, 1}, 
    {0, 0, 1, 1}, 
    {0, 0, 0, 1}, 
}};

void dJydy_colptrs_model_robertson_py(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_model_robertson_py_[index]));
}
} // namespace model_model_robertson_py
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_robertson_py {

static constexpr std::array<std::array<sunindextype, 1>, 3> dJydy_rowvals_model_robertson_py_ = {{
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_model_robertson_py(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_model_robertson_py_[index]));
}
} // namespace model_model_robertson_py
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include <sundials/sundials_types.h>
#include <gsl/gsl-lite.hpp>

namespace amici {
namespace model_model_robertson_py {

void dJydy_model_robertson_py(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    const realtype obs_x1_ = y[0];
    const realtype obs_x2_ = y[1];
    const realtype obs_x3_ = y[2];
    const realtype sigma_obs_x1_ = sigmay[0];
    const realtype sigma_obs_x2_ = sigmay[1];
    const realtype sigma_obs_x3_ = sigmay[2];
    const realtype mobs_x1_ = my[0];
    const realtype mobs_x2_ = my[1];
    const realtype mobs_x3_ = my[2];

    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mobs_x1_ + 1.0*obs_x1_)/std::pow(sigma_obs_x1_, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*mobs_x2_ + 1.0*obs_x2_)/std::pow(sigma_obs_x2_, 2);
            break;
        case 2:
            dJydy[0] = (-1.0*mobs_x3_ + 1.0*obs_x3_)/std::pow(sigma_obs_x3_, 2);
            break;
    }
}

} // namespace model_model_robertson_py
} // namespace amici
