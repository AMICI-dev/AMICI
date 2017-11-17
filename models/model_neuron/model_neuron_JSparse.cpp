
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void JSparse_model_neuron(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  J[0] = x[0]*(2.0/2.5E1)+5.0;
  J[1] = p[0]*p[1];
  J[2] = -1.0;
  J[3] = -p[0];
}

