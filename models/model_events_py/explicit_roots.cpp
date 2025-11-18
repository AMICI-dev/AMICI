#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include <vector>
#include "p.h"
#include "k.h"

namespace amici {
namespace model_model_events_py {

std::vector<std::vector<realtype>> explicit_roots_model_events_py(const realtype *p, const realtype *k, const realtype *w){
    return {
        {p4},
        {p4},
        {4},
        {4}
    };
}

} // namespace model_model_events_py
} // namespace amici
