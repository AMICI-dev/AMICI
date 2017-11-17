
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void Jv_model_dirac(realtype *Jv, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *v, const realtype *w, const realtype *dwdx) {
  Jv[0] = -p[0]*v[0];
  Jv[1] = p[2]*v[0]-p[3]*v[1];
}

