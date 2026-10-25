#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_jakstat_adjoint_py {

static constexpr std::array<std::array<sunindextype, 4>, 3> dJydy_colptrs_model_jakstat_adjoint_py_ = {{
    {0, 1, 1, 1}, 
    {0, 0, 1, 1}, 
    {0, 0, 0, 1}, 
}};

void dJydy_colptrs_model_jakstat_adjoint_py(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_model_jakstat_adjoint_py_[index]));
}
} // namespace model_model_jakstat_adjoint_py
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_jakstat_adjoint_py {

static constexpr std::array<std::array<sunindextype, 1>, 3> dJydy_rowvals_model_jakstat_adjoint_py_ = {{
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_model_jakstat_adjoint_py(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_model_jakstat_adjoint_py_[index]));
}
} // namespace model_model_jakstat_adjoint_py
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include <sundials/sundials_types.h>
#include <gsl/gsl-lite.hpp>

namespace amici {
namespace model_model_jakstat_adjoint_py {

void dJydy_model_jakstat_adjoint_py(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    const realtype obs_pSTAT_ = y[0];
    const realtype obs_tSTAT_ = y[1];
    const realtype obs_spline_ = y[2];
    const realtype sigma_obs_pSTAT_ = sigmay[0];
    const realtype sigma_obs_tSTAT_ = sigmay[1];
    const realtype sigma_obs_spline_ = sigmay[2];
    const realtype mobs_pSTAT_ = my[0];
    const realtype mobs_tSTAT_ = my[1];
    const realtype mobs_spline_ = my[2];

    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mobs_pSTAT_ + 1.0*obs_pSTAT_)/std::pow(sigma_obs_pSTAT_, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*mobs_tSTAT_ + 1.0*obs_tSTAT_)/std::pow(sigma_obs_tSTAT_, 2);
            break;
        case 2:
            dJydy[0] = (-1.0*mobs_spline_ + 1.0*obs_spline_)/std::pow(sigma_obs_spline_, 2);
            break;
    }
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
