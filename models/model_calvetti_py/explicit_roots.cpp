#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include <vector>
#include "k.h"

namespace amici {
namespace model_model_calvetti_py {

std::vector<std::vector<realtype>> explicit_roots_model_calvetti_py(const realtype *p, const realtype *k){
    return {
        {10},
        {10},
        {12},
        {12}
    };
}

} // namespace model_model_calvetti_py
} // namespace amici
