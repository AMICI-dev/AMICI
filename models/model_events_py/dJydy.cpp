#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_events_py {

static constexpr std::array<std::array<sunindextype, 2>, 1> dJydy_colptrs_model_events_py_ = {{
    {0, 1}, 
}};

void dJydy_colptrs_model_events_py(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_model_events_py_[index]));
}
} // namespace model_model_events_py
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_model_events_py {

static constexpr std::array<std::array<sunindextype, 1>, 1> dJydy_rowvals_model_events_py_ = {{
    {0}, 
}};

void dJydy_rowvals_model_events_py(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_model_events_py_[index]));
}
} // namespace model_model_events_py
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"
#include "sigmay.h"
#include "my.h"
#include "dJydy.h"

namespace amici {
namespace model_model_events_py {

void dJydy_model_events_py(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*my1 + 1.0*y1)/std::pow(sigma_y1, 2);
            break;
    }
}

} // namespace model_model_events_py
} // namespace amici
