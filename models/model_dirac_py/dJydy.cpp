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
#include "p.h"
#include "y.h"
#include "sigmay.h"
#include "my.h"
#include "dJydy.h"

namespace amici {
namespace model_model_dirac_py {

void dJydy_model_dirac_py(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mobs_x2 + 1.0*obs_x2)/std::pow(sigma_obs_x2, 2);
            break;
    }
}

} // namespace model_model_dirac_py
} // namespace amici
