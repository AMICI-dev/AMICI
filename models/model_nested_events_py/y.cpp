#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>

namespace amici {
namespace model_model_nested_events_py {

void y_model_nested_events_py(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    const realtype Virus_ = x[0];

    y[0] = Virus_;
}

} // namespace model_model_nested_events_py
} // namespace amici
