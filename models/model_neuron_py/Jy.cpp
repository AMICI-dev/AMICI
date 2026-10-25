#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_neuron_py {

void Jy_model_neuron_py(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    const realtype y1_ = y[0];
    const realtype sigma_y1_ = sigmay[0];
    const realtype my1_ = my[0];

    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_y1_, 2)) + 0.5*std::pow(-my1_ + y1_, 2)/std::pow(sigma_y1_, 2);
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
