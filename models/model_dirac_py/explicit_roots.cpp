#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include <vector>

namespace amici {
namespace model_model_dirac_py {

std::vector<std::vector<realtype>> explicit_roots_model_dirac_py(const realtype *x, const realtype *p, const realtype *k, const realtype *w){
    const realtype p2_ = p[1];

    return {
        {p2_}
    };
}

} // namespace model_model_dirac_py
} // namespace amici
