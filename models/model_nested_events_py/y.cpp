#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "x.h"
#include "p.h"
#include "h.h"

namespace amici {
namespace model_model_nested_events_py {

void y_model_nested_events_py(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = Virus;
}

} // namespace model_model_nested_events_py
} // namespace amici
