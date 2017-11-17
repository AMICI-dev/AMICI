
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void JSparseB_model_neuron(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  JB[0] = x[0]*(-2.0/2.5E1)-5.0;
  JB[1] = 1.0;
  JB[2] = -p[0]*p[1];
  JB[3] = p[0];
}

