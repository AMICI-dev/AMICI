#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include <gsl/gsl-lite.hpp>

namespace amici {
namespace model_model_events_py {

void x0_fixedParameters_model_events_py(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs){
    const realtype k1_ = k[0];
    const realtype k2_ = k[1];
    const realtype k3_ = k[2];

    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 0) != reinitialization_state_idxs.cend())
        x0_fixedParameters[0] = k1_;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 1) != reinitialization_state_idxs.cend())
        x0_fixedParameters[1] = k2_;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 2) != reinitialization_state_idxs.cend())
        x0_fixedParameters[2] = k3_;
}

} // namespace model_model_events_py
} // namespace amici
