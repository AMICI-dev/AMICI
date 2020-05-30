
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
#include <sunmatrix/sunmatrix_sparse.h> //SUNMatrixContent_Sparse definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_steadystate{

void JSparseB_model_steadystate(SUNMatrixContent_Sparse JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  JSparseB->indexvals[0] = 0;
  JSparseB->indexvals[1] = 1;
  JSparseB->indexvals[2] = 2;
  JSparseB->indexvals[3] = 0;
  JSparseB->indexvals[4] = 1;
  JSparseB->indexvals[5] = 2;
  JSparseB->indexvals[6] = 0;
  JSparseB->indexvals[7] = 1;
  JSparseB->indexvals[8] = 2;
  JSparseB->indexptrs[0] = 0;
  JSparseB->indexptrs[1] = 3;
  JSparseB->indexptrs[2] = 6;
  JSparseB->indexptrs[3] = 9;
  JSparseB->data[0] = p[1]*x[1]+p[0]*dwdx[0]*2.0;
  JSparseB->data[1] = p[2]*-2.0+p[1]*x[0];
  JSparseB->data[2] = -dwdx[1];
  JSparseB->data[3] = p[1]*x[1]-p[0]*dwdx[0];
  JSparseB->data[4] = p[2]+p[1]*x[0];
  JSparseB->data[5] = -dwdx[1];
  JSparseB->data[6] = -p[1]*x[1];
  JSparseB->data[7] = -p[1]*x[0];
  JSparseB->data[8] = k[3]+dwdx[1];
}

} // namespace model_model_steadystate

} // namespace amici

