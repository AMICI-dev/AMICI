
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
#include <sunmatrix/sunmatrix_sparse.h> //SUNMatrixContent_Sparse definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron_o2{

void JSparse_model_neuron_o2(SUNMatrixContent_Sparse JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  JSparse->indexvals[0] = 0;
  JSparse->indexvals[1] = 1;
  JSparse->indexvals[2] = 2;
  JSparse->indexvals[3] = 3;
  JSparse->indexvals[4] = 4;
  JSparse->indexvals[5] = 5;
  JSparse->indexvals[6] = 6;
  JSparse->indexvals[7] = 8;
  JSparse->indexvals[8] = 0;
  JSparse->indexvals[9] = 1;
  JSparse->indexvals[10] = 3;
  JSparse->indexvals[11] = 2;
  JSparse->indexvals[12] = 3;
  JSparse->indexvals[13] = 2;
  JSparse->indexvals[14] = 3;
  JSparse->indexvals[15] = 4;
  JSparse->indexvals[16] = 5;
  JSparse->indexvals[17] = 4;
  JSparse->indexvals[18] = 5;
  JSparse->indexvals[19] = 6;
  JSparse->indexvals[20] = 7;
  JSparse->indexvals[21] = 6;
  JSparse->indexvals[22] = 7;
  JSparse->indexvals[23] = 8;
  JSparse->indexvals[24] = 9;
  JSparse->indexvals[25] = 8;
  JSparse->indexvals[26] = 9;
  JSparse->indexptrs[0] = 0;
  JSparse->indexptrs[1] = 8;
  JSparse->indexptrs[2] = 11;
  JSparse->indexptrs[3] = 13;
  JSparse->indexptrs[4] = 15;
  JSparse->indexptrs[5] = 17;
  JSparse->indexptrs[6] = 19;
  JSparse->indexptrs[7] = 21;
  JSparse->indexptrs[8] = 23;
  JSparse->indexptrs[9] = 25;
  JSparse->indexptrs[10] = 27;
  JSparse->data[0] = x[0]*(2.0/2.5E1)+5.0;
  JSparse->data[1] = p[0]*p[1];
  JSparse->data[2] = x[2]*dwdx[1];
  JSparse->data[3] = p[1];
  JSparse->data[4] = x[4]*dwdx[1];
  JSparse->data[5] = p[0];
  JSparse->data[6] = x[6]*dwdx[1];
  JSparse->data[7] = x[8]*dwdx[1];
  JSparse->data[8] = -1.0;
  JSparse->data[9] = -p[0];
  JSparse->data[10] = -1.0;
  JSparse->data[11] = w[1];
  JSparse->data[12] = p[0]*p[1];
  JSparse->data[13] = -1.0;
  JSparse->data[14] = -p[0];
  JSparse->data[15] = w[1];
  JSparse->data[16] = p[0]*p[1];
  JSparse->data[17] = -1.0;
  JSparse->data[18] = -p[0];
  JSparse->data[19] = w[1];
  JSparse->data[20] = p[0]*p[1];
  JSparse->data[21] = -1.0;
  JSparse->data[22] = -p[0];
  JSparse->data[23] = w[1];
  JSparse->data[24] = p[0]*p[1];
  JSparse->data[25] = -1.0;
  JSparse->data[26] = -p[0];
}

} // namespace model_model_neuron_o2

} // namespace amici

