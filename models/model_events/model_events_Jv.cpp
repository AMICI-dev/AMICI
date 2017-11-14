
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void Jv_model_events(realtype *Jv, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *v, const realtype *w, const realtype *dwdx) {
  Jv[0] = -h[3]*p[0]*v[0];
  Jv[1] = -p[2]*v[1]+p[1]*v[0]*exp(t*(-1.0/1.0E1));
  Jv[2] = -v[2];
}

