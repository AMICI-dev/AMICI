#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include <vector>

namespace amici {
namespace model_model_events_py {

std::vector<std::vector<realtype>> explicit_roots_model_events_py(const realtype *x, const realtype *p, const realtype *k, const realtype *w){
    const realtype p4_ = p[3];

    return {
        {p4_},
        {p4_},
        {4},
        {4}
    };
}

} // namespace model_model_events_py
} // namespace amici
