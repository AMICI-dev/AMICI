
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void x0_model_jakstat_adjoint_o2(realtype *x0, const realtype t, const realtype *p, const realtype *k) {
  x0[0] = p[4];
  x0[45] = 1.0;
}

