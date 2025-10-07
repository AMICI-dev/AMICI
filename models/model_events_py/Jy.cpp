#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "p.h"
#include "k.h"
#include "y.h"
#include "sigmay.h"
#include "my.h"

namespace amici {
namespace model_model_events_py {

void Jy_model_events_py(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_y1, 2)) + 0.5*std::pow(-my1 + y1, 2)/std::pow(sigma_y1, 2);
            break;
    }
}

} // namespace model_model_events_py
} // namespace amici
