
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
#include <sunmatrix/sunmatrix_sparse.h> //SUNMatrixContent_Sparse definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_calvetti{

void JSparse_model_calvetti(SUNMatrixContent_Sparse JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *dx, const realtype *w, const realtype *dwdx) {
  JSparse->indexvals[0] = 0;
  JSparse->indexvals[1] = 3;
  JSparse->indexvals[2] = 0;
  JSparse->indexvals[3] = 1;
  JSparse->indexvals[4] = 3;
  JSparse->indexvals[5] = 4;
  JSparse->indexvals[6] = 0;
  JSparse->indexvals[7] = 1;
  JSparse->indexvals[8] = 2;
  JSparse->indexvals[9] = 3;
  JSparse->indexvals[10] = 4;
  JSparse->indexvals[11] = 5;
  JSparse->indexvals[12] = 0;
  JSparse->indexvals[13] = 3;
  JSparse->indexvals[14] = 4;
  JSparse->indexvals[15] = 0;
  JSparse->indexvals[16] = 1;
  JSparse->indexvals[17] = 3;
  JSparse->indexvals[18] = 4;
  JSparse->indexvals[19] = 5;
  JSparse->indexvals[20] = 0;
  JSparse->indexvals[21] = 1;
  JSparse->indexvals[22] = 2;
  JSparse->indexvals[23] = 3;
  JSparse->indexvals[24] = 4;
  JSparse->indexvals[25] = 5;
  JSparse->indexptrs[0] = 0;
  JSparse->indexptrs[1] = 2;
  JSparse->indexptrs[2] = 6;
  JSparse->indexptrs[3] = 12;
  JSparse->indexptrs[4] = 15;
  JSparse->indexptrs[5] = 20;
  JSparse->indexptrs[6] = 26;
  JSparse->data[0] = -cj-w[0]*(1.0E2/8.99E2)+dwdx[6];
  JSparse->data[1] = w[0]*(1.0E2/8.99E2)-dwdx[6]+dwdx[7]*(w[23]*w[24]*2.0-w[20]*w[23]*w[24])-w[23]*w[24]*w[25]*dwdx[4];
  JSparse->data[2] = dwdx[16];
  JSparse->data[3] = -cj-w[6]*1.202935161794779E-2+dwdx[19];
  JSparse->data[4] = -dwdx[16]-w[23]*w[24]*w[25]*dwdx[14];
  JSparse->data[5] = w[6]*1.202935161794779E-2-dwdx[19];
  JSparse->data[6] = dwdx[28];
  JSparse->data[7] = dwdx[31];
  JSparse->data[8] = -cj-w[31]*8.196729508204918E-9+dwdx[32];
  JSparse->data[9] = -dwdx[28]-w[23]*w[24]*w[25]*dwdx[26];
  JSparse->data[10] = -dwdx[31];
  JSparse->data[11] = w[31]*8.196729508204918E-9-dwdx[32];
  JSparse->data[12] = dwdx[36];
  JSparse->data[13] = -dwdx[36]-w[23]*w[24]*w[25]*dwdx[34]-1.0;
  JSparse->data[14] = 1.0;
  JSparse->data[15] = dwdx[40];
  JSparse->data[16] = dwdx[43];
  JSparse->data[17] = -dwdx[40]-w[23]*w[24]*w[25]*dwdx[38];
  JSparse->data[18] = -dwdx[43]-1.0;
  JSparse->data[19] = 1.0;
  JSparse->data[20] = dwdx[47];
  JSparse->data[21] = dwdx[50];
  JSparse->data[22] = dwdx[52];
  JSparse->data[23] = -dwdx[47]-w[23]*w[24]*w[25]*dwdx[45];
  JSparse->data[24] = -dwdx[50];
  JSparse->data[25] = -dwdx[52]-1.0;
}

} // namespace model_model_calvetti

} // namespace amici

