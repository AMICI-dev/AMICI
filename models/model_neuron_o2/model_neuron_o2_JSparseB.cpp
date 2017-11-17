
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void JSparseB_model_neuron_o2(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  JB[0] = x[0]*(-2.0/2.5E1)-5.0;
  JB[1] = 1.0;
  JB[2] = -p[0]*p[1];
  JB[3] = p[0];
  JB[5] = -w[1];
  JB[6] = 1.0;
  JB[9] = -p[0]*p[1];
  JB[10] = p[0];
  JB[12] = -w[1];
  JB[13] = 1.0;
  JB[15] = -p[0]*p[1];
  JB[16] = p[0];
  JB[18] = -w[1];
  JB[19] = 1.0;
  JB[20] = -p[0]*p[1];
  JB[21] = p[0];
  JB[23] = -w[1];
  JB[24] = 1.0;
  JB[25] = -p[0]*p[1];
  JB[26] = p[0];
}

