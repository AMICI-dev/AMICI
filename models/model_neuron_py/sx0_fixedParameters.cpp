#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include <gsl/gsl-lite.hpp>
#include "p.h"
#include "k.h"

namespace amici {
namespace model_model_neuron_py {

void sx0_fixedParameters_model_neuron_py(realtype *sx0_fixedParameters, const realtype t, const realtype *x0, const realtype *p, const realtype *k, const int ip, gsl::span<const int> reinitialization_state_idxs){
    static const std::array<int, 2> _x0_fixedParameters_idxs = {
        0, 1
    };
    for(auto idx: reinitialization_state_idxs) {
        if(std::find(_x0_fixedParameters_idxs.cbegin(), _x0_fixedParameters_idxs.cend(), idx) != _x0_fixedParameters_idxs.cend())
            sx0_fixedParameters[idx] = 0.0;
    }
    switch(ip) {
        case 1:
            if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 1) != reinitialization_state_idxs.cend())
                sx0_fixedParameters[1] = v0;
            break;
    }
}

} // namespace model_model_neuron_py
} // namespace amici
