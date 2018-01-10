
#include <include/symbolic_functions.h>
#include <include/amici_defines.h> //realtype definition
#include <sundials/sundials_sparse.h> //SlsMat definition
typedef amici::realtype realtype;
#include <cmath> 

void JSparseB_model_neuron_o2(SlsMat JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  JSparseB->indexvals[0] = 0;
  JSparseB->indexvals[1] = 1;
  JSparseB->indexvals[2] = 0;
  JSparseB->indexvals[3] = 1;
  JSparseB->indexvals[4] = 0;
  JSparseB->indexvals[5] = 2;
  JSparseB->indexvals[6] = 3;
  JSparseB->indexvals[7] = 0;
  JSparseB->indexvals[8] = 1;
  JSparseB->indexvals[9] = 2;
  JSparseB->indexvals[10] = 3;
  JSparseB->indexvals[11] = 0;
  JSparseB->indexvals[12] = 4;
  JSparseB->indexvals[13] = 5;
  JSparseB->indexvals[14] = 0;
  JSparseB->indexvals[15] = 4;
  JSparseB->indexvals[16] = 5;
  JSparseB->indexvals[17] = 0;
  JSparseB->indexvals[18] = 6;
  JSparseB->indexvals[19] = 7;
  JSparseB->indexvals[20] = 6;
  JSparseB->indexvals[21] = 7;
  JSparseB->indexvals[22] = 0;
  JSparseB->indexvals[23] = 8;
  JSparseB->indexvals[24] = 9;
  JSparseB->indexvals[25] = 8;
  JSparseB->indexvals[26] = 9;
  JSparseB->indexptrs[0] = 0;
  JSparseB->indexptrs[1] = 2;
  JSparseB->indexptrs[2] = 4;
  JSparseB->indexptrs[3] = 7;
  JSparseB->indexptrs[4] = 11;
  JSparseB->indexptrs[5] = 14;
  JSparseB->indexptrs[6] = 17;
  JSparseB->indexptrs[7] = 20;
  JSparseB->indexptrs[8] = 22;
  JSparseB->indexptrs[9] = 25;
  JSparseB->indexptrs[10] = 27;
  JSparseB->data[0] = x[0]*(-2.0/2.5E1)-5.0;
  JSparseB->data[1] = 1.0;
  JSparseB->data[2] = -p[0]*p[1];
  JSparseB->data[3] = p[0];
  JSparseB->data[5] = -w[1];
  JSparseB->data[6] = 1.0;
  JSparseB->data[9] = -p[0]*p[1];
  JSparseB->data[10] = p[0];
  JSparseB->data[12] = -w[1];
  JSparseB->data[13] = 1.0;
  JSparseB->data[15] = -p[0]*p[1];
  JSparseB->data[16] = p[0];
  JSparseB->data[18] = -w[1];
  JSparseB->data[19] = 1.0;
  JSparseB->data[20] = -p[0]*p[1];
  JSparseB->data[21] = p[0];
  JSparseB->data[23] = -w[1];
  JSparseB->data[24] = 1.0;
  JSparseB->data[25] = -p[0]*p[1];
  JSparseB->data[26] = p[0];
}

