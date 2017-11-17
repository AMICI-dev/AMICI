
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void JSparse_model_neuron_o2(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  J[0] = x[0]*(2.0/2.5E1)+5.0;
  J[1] = p[0]*p[1];
  J[2] = x[2]*dwdx[1];
  J[3] = p[1];
  J[4] = x[4]*dwdx[1];
  J[5] = p[0];
  J[6] = x[6]*dwdx[1];
  J[7] = x[8]*dwdx[1];
  J[8] = -1.0;
  J[9] = -p[0];
  J[10] = -1.0;
  J[11] = w[1];
  J[12] = p[0]*p[1];
  J[13] = -1.0;
  J[14] = -p[0];
  J[15] = w[1];
  J[16] = p[0]*p[1];
  J[17] = -1.0;
  J[18] = -p[0];
  J[19] = w[1];
  J[20] = p[0]*p[1];
  J[21] = -1.0;
  J[22] = -p[0];
  J[23] = w[1];
  J[24] = p[0]*p[1];
  J[25] = -1.0;
  J[26] = -p[0];
}

