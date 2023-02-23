
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
#include <sunmatrix/sunmatrix_sparse.h> //SUNMatrixContent_Sparse definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_steadystate{

void JSparse_model_steadystate(SUNMatrixContent_Sparse JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) {
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
  JSparse->data[0] = -p[1]*x[1]-p[0]*dwdx[0]*2.0;
  JSparse->data[1] = -p[1]*x[1]+p[0]*dwdx[0];
  JSparse->data[2] = p[1]*x[1];
  JSparse->data[3] = p[2]*2.0-p[1]*x[0];
  JSparse->data[4] = -p[2]-p[1]*x[0];
  JSparse->data[5] = p[1]*x[0];
  JSparse->data[6] = dwdx[1];
  JSparse->data[7] = dwdx[1];
  JSparse->data[8] = -k[3]-dwdx[1];
}

} // namespace model_model_steadystate

} // namespace amici

