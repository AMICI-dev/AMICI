#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include <gsl/gsl-lite.hpp>

namespace amici {
namespace model_model_neuron_py {

void x0_fixedParameters_model_neuron_py(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs){
    const realtype b_ = p[1];
    const realtype v0_ = k[0];

    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 0) != reinitialization_state_idxs.cend())
        x0_fixedParameters[0] = v0_;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 1) != reinitialization_state_idxs.cend())
        x0_fixedParameters[1] = b_*v0_;
}

} // namespace model_model_neuron_py
} // namespace amici
