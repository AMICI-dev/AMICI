
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void Jv_model_neuron_o2(realtype *Jv, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *v, const realtype *w, const realtype *dwdx) {
  Jv[0] = -v[1]+v[0]*(x[0]*(2.0/2.5E1)+5.0);
  Jv[1] = -p[0]*v[1]+p[0]*p[1]*v[0];
  Jv[2] = -v[3]+v[2]*w[1]+v[0]*x[2]*dwdx[1];
  Jv[3] = -v[1]+p[1]*v[0]-p[0]*v[3]+p[0]*p[1]*v[2];
  Jv[4] = -v[5]+w[1]*v[4]+v[0]*x[4]*dwdx[1];
  Jv[5] = p[0]*v[0]-p[0]*v[5]+p[0]*p[1]*v[4];
  Jv[6] = -v[7]+w[1]*v[6]+v[0]*x[6]*dwdx[1];
  Jv[7] = -p[0]*v[7]+p[0]*p[1]*v[6];
  Jv[8] = -v[9]+w[1]*v[8]+v[0]*x[8]*dwdx[1];
  Jv[9] = -p[0]*v[9]+p[0]*p[1]*v[8];
}

