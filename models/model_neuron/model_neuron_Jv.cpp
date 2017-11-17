
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void Jv_model_neuron(realtype *Jv, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *v, const realtype *w, const realtype *dwdx) {
  Jv[0] = -v[1]+v[0]*(x[0]*(2.0/2.5E1)+5.0);
  Jv[1] = -p[0]*v[1]+p[0]*p[1]*v[0];
}

