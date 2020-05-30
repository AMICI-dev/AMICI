
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
#include <sunmatrix/sunmatrix_sparse.h> //SUNMatrixContent_Sparse definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_jakstat_adjoint{

void JSparse_model_jakstat_adjoint(SUNMatrixContent_Sparse JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  JSparse->indexvals[0] = 0;
  JSparse->indexvals[1] = 1;
  JSparse->indexvals[2] = 1;
  JSparse->indexvals[3] = 2;
  JSparse->indexvals[4] = 2;
  JSparse->indexvals[5] = 3;
  JSparse->indexvals[6] = 3;
  JSparse->indexvals[7] = 4;
  JSparse->indexvals[8] = 4;
  JSparse->indexvals[9] = 5;
  JSparse->indexvals[10] = 5;
  JSparse->indexvals[11] = 6;
  JSparse->indexvals[12] = 6;
  JSparse->indexvals[13] = 7;
  JSparse->indexvals[14] = 7;
  JSparse->indexvals[15] = 8;
  JSparse->indexvals[16] = 0;
  JSparse->indexvals[17] = 8;
  JSparse->indexptrs[0] = 0;
  JSparse->indexptrs[1] = 2;
  JSparse->indexptrs[2] = 4;
  JSparse->indexptrs[3] = 6;
  JSparse->indexptrs[4] = 8;
  JSparse->indexptrs[5] = 10;
  JSparse->indexptrs[6] = 12;
  JSparse->indexptrs[7] = 14;
  JSparse->indexptrs[8] = 16;
  JSparse->indexptrs[9] = 18;
  JSparse->data[0] = -p[0]*w[0];
  JSparse->data[1] = p[0]*w[0];
  JSparse->data[2] = p[1]*dwdx[0]*-2.0;
  JSparse->data[3] = p[1]*dwdx[0];
  JSparse->data[4] = -p[2];
  JSparse->data[5] = (k[0]*p[2])/k[1];
  JSparse->data[6] = -p[3];
  JSparse->data[7] = p[3]*2.0;
  JSparse->data[8] = -p[3];
  JSparse->data[9] = p[3];
  JSparse->data[10] = -p[3];
  JSparse->data[11] = p[3];
  JSparse->data[12] = -p[3];
  JSparse->data[13] = p[3];
  JSparse->data[14] = -p[3];
  JSparse->data[15] = p[3];
  JSparse->data[16] = (k[1]*p[3])/k[0];
  JSparse->data[17] = -p[3];
}

} // namespace model_model_jakstat_adjoint

} // namespace amici

