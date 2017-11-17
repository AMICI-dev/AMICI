
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <sundials/sundials_sparse.h> //SlsMat definition
#include <cmath> 

void JSparse_model_robertson(SlsMat JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *dx, const realtype *w, const realtype *dwdx) {
  JSparse->indexvals[0] = 0;
  JSparse->indexvals[1] = 1;
  JSparse->indexvals[2] = 2;
  JSparse->indexvals[3] = 0;
  JSparse->indexvals[4] = 1;
  JSparse->indexvals[5] = 2;
  JSparse->indexvals[6] = 0;
  JSparse->indexvals[7] = 1;
  JSparse->indexvals[8] = 2;
  JSparse->indexptrs[0] = 0;
  JSparse->indexptrs[1] = 3;
  JSparse->indexptrs[2] = 6;
  JSparse->indexptrs[3] = 9;
  JSparse->data[0] = -cj-p[0];
  JSparse->data[1] = p[0];
  JSparse->data[2] = 1.0;
  JSparse->data[3] = dwdx[0];
  JSparse->data[4] = -cj-dwdx[0]-p[2]*x[1]*2.0;
  JSparse->data[5] = 1.0;
  JSparse->data[6] = dwdx[1];
  JSparse->data[7] = -dwdx[1];
  JSparse->data[8] = 1.0;
}

