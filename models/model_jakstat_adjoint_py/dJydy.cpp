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
#include "p.h"
#include "k.h"
#include "y.h"
#include "sigmay.h"
#include "my.h"
#include "dJydy.h"

namespace amici {
namespace model_model_jakstat_adjoint_py {

void dJydy_model_jakstat_adjoint_py(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mobs_pSTAT + 1.0*obs_pSTAT)/std::pow(sigma_obs_pSTAT, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*mobs_tSTAT + 1.0*obs_tSTAT)/std::pow(sigma_obs_tSTAT, 2);
            break;
        case 2:
            dJydy[0] = (-1.0*mobs_spline + 1.0*obs_spline)/std::pow(sigma_obs_spline, 2);
            break;
    }
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
