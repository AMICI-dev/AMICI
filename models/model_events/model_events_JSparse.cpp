
#include <include/symbolic_functions.h>
#include <include/amici_defines.h> //realtype definition
#include <sundials/sundials_sparse.h> //SlsMat definition
typedef amici::realtype realtype;
#include <cmath> 

void JSparse_model_events(SlsMat JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  JSparse->indexvals[0] = 0;
  JSparse->indexvals[1] = 1;
  JSparse->indexvals[2] = 1;
  JSparse->indexvals[3] = 2;
  JSparse->indexptrs[0] = 0;
  JSparse->indexptrs[1] = 2;
  JSparse->indexptrs[2] = 3;
  JSparse->indexptrs[3] = 4;
  JSparse->data[0] = -h[3]*p[0];
  JSparse->data[1] = p[1]*exp(t*(-1.0/1.0E1));
  JSparse->data[2] = -p[2];
  JSparse->data[3] = -1.0;
}

