#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_events_py {

void dJydsigmay_model_events_py(realtype *dJydsigmay, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    const realtype y1_ = y[0];
    const realtype sigma_y1_ = sigmay[0];
    const realtype my1_ = my[0];

    switch(iy) {
        case 0:
            dJydsigmay[0] = 1.0/sigma_y1_ - 1.0*std::pow(-my1_ + y1_, 2)/std::pow(sigma_y1_, 3);
            break;
    }
}

} // namespace model_model_events_py
} // namespace amici
