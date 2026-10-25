#include "amici/symbolic_functions.h"
#include "amici/defines.h"

#include <algorithm>
#include "amici/splinefunctions.h"
#include <vector>

namespace amici {
namespace model_model_jakstat_adjoint_py {

std::vector<HermiteSpline> create_splines_model_jakstat_adjoint_py(const realtype *p, const realtype *k){
    const realtype sp1_ = p[5];
    const realtype sp2_ = p[6];
    const realtype sp3_ = p[7];
    const realtype sp4_ = p[8];
    const realtype sp5_ = p[9];

    return {
        HermiteSpline(
            {0, 5, 10, 20, 60}, 
            {sp1_, sp2_, sp3_, sp4_, sp5_}, 
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
