#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "amici/splinefunctions.h"
#include <vector>
#include "p.h"
#include "k.h"

namespace amici {
namespace model_model_events_py {

std::vector<HermiteSpline> create_splines_model_events_py(const realtype *p, const realtype *k){
    return {};
}

} // namespace model_model_events_py
} // namespace amici
