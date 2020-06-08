
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
#include <sunmatrix/sunmatrix_sparse.h> //SUNMatrixContent_Sparse definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_jakstat_adjoint{

void JSparseB_model_jakstat_adjoint(SUNMatrixContent_Sparse JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  JSparseB->indexvals[0] = 0;
  JSparseB->indexvals[1] = 8;
  JSparseB->indexvals[2] = 0;
  JSparseB->indexvals[3] = 1;
  JSparseB->indexvals[4] = 1;
  JSparseB->indexvals[5] = 2;
  JSparseB->indexvals[6] = 2;
  JSparseB->indexvals[7] = 3;
  JSparseB->indexvals[8] = 3;
  JSparseB->indexvals[9] = 4;
  JSparseB->indexvals[10] = 4;
  JSparseB->indexvals[11] = 5;
  JSparseB->indexvals[12] = 5;
  JSparseB->indexvals[13] = 6;
  JSparseB->indexvals[14] = 6;
  JSparseB->indexvals[15] = 7;
  JSparseB->indexvals[16] = 7;
  JSparseB->indexvals[17] = 8;
  JSparseB->indexptrs[0] = 0;
  JSparseB->indexptrs[1] = 2;
  JSparseB->indexptrs[2] = 4;
  JSparseB->indexptrs[3] = 6;
  JSparseB->indexptrs[4] = 8;
  JSparseB->indexptrs[5] = 10;
  JSparseB->indexptrs[6] = 12;
  JSparseB->indexptrs[7] = 14;
  JSparseB->indexptrs[8] = 16;
  JSparseB->indexptrs[9] = 18;
  JSparseB->data[0] = p[0]*w[0];
  JSparseB->data[1] = -(k[1]*p[3])/k[0];
  JSparseB->data[2] = -p[0]*w[0];
  JSparseB->data[3] = p[1]*dwdx[0]*2.0;
  JSparseB->data[4] = -p[1]*dwdx[0];
  JSparseB->data[5] = p[2];
  JSparseB->data[6] = -(k[0]*p[2])/k[1];
  JSparseB->data[7] = p[3];
  JSparseB->data[8] = p[3]*-2.0;
  JSparseB->data[9] = p[3];
  JSparseB->data[10] = -p[3];
  JSparseB->data[11] = p[3];
  JSparseB->data[12] = -p[3];
  JSparseB->data[13] = p[3];
  JSparseB->data[14] = -p[3];
  JSparseB->data[15] = p[3];
  JSparseB->data[16] = -p[3];
  JSparseB->data[17] = p[3];
}

} // namespace model_model_jakstat_adjoint

} // namespace amici

