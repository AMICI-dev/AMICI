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
namespace model_model_jakstat_adjoint_py {

std::vector<HermiteSpline> create_splines_model_jakstat_adjoint_py(const realtype *p, const realtype *k){
    return {
        HermiteSpline(
            {0, 5, 10, 20, 60}, 
            {sp1, sp2, sp3, sp4, sp5}, 
            {},
            SplineBoundaryCondition::zeroDerivative, 
            SplineBoundaryCondition::zeroDerivative, 
            SplineExtrapolation::constant, 
            SplineExtrapolation::constant, 
            true, false, true
        ),
    };
}

} // namespace model_model_jakstat_adjoint_py
} // namespace amici
