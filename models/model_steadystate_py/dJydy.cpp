#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_steadystate_py {

static constexpr std::array<std::array<sunindextype, 4>, 3> dJydy_colptrs_model_steadystate_py_ = {{
    {0, 1, 1, 1}, 
    {0, 0, 1, 1}, 
    {0, 0, 0, 1}, 
}};

void dJydy_colptrs_model_steadystate_py(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_model_steadystate_py_[index]));
}
} // namespace model_model_steadystate_py
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_steadystate_py {

static constexpr std::array<std::array<sunindextype, 1>, 3> dJydy_rowvals_model_steadystate_py_ = {{
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_model_steadystate_py(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_model_steadystate_py_[index]));
}
} // namespace model_model_steadystate_py
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include <sundials/sundials_types.h>
#include <gsl/gsl-lite.hpp>
#include "p.h"
#include "k.h"
#include "y.h"
#include "sigmay.h"
#include "my.h"
#include "dJydy.h"

namespace amici {
namespace model_model_steadystate_py {

void dJydy_model_steadystate_py(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mobs_x1 + 1.0*obs_x1)/std::pow(sigma_obs_x1, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*mobs_x2 + 1.0*obs_x2)/std::pow(sigma_obs_x2, 2);
            break;
        case 2:
            dJydy[0] = (-1.0*mobs_x3 + 1.0*obs_x3)/std::pow(sigma_obs_x3, 2);
            break;
    }
}

} // namespace model_model_steadystate_py
} // namespace amici
