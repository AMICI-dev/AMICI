#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_dirac_py {

static constexpr std::array<std::array<sunindextype, 2>, 1> dJydy_colptrs_model_dirac_py_ = {{
    {0, 1}, 
}};

void dJydy_colptrs_model_dirac_py(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_model_dirac_py_[index]));
}
} // namespace model_model_dirac_py
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_dirac_py {

static constexpr std::array<std::array<sunindextype, 1>, 1> dJydy_rowvals_model_dirac_py_ = {{
    {0}, 
}};

void dJydy_rowvals_model_dirac_py(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_model_dirac_py_[index]));
}
} // namespace model_model_dirac_py
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include <sundials/sundials_types.h>
#include <gsl/gsl-lite.hpp>

namespace amici {
namespace model_model_dirac_py {

void dJydy_model_dirac_py(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    const realtype obs_x2_ = y[0];
    const realtype sigma_obs_x2_ = sigmay[0];
    const realtype mobs_x2_ = my[0];

    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mobs_x2_ + 1.0*obs_x2_)/std::pow(sigma_obs_x2_, 2);
            break;
    }
}

} // namespace model_model_dirac_py
} // namespace amici
