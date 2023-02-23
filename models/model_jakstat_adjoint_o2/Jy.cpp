
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_jakstat_adjoint_o2{

void Jy_model_jakstat_adjoint_o2(double *nllh, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) {
switch(iy){
    case 0:
  nllh[0] = amici::log((sigmay[0]*sigmay[0])*3.141592653589793*2.0)*5.0E-1+1.0/(sigmay[0]*sigmay[0])*pow(my[0]-y[0],2.0)*5.0E-1;
  nllh[1] = y[3]*1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  nllh[2] = y[6]*1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  nllh[3] = y[9]*1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  nllh[4] = y[12]*1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  nllh[5] = y[15]*1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  nllh[6] = y[18]*1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  nllh[7] = y[21]*1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  nllh[8] = y[24]*1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  nllh[9] = y[27]*1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  nllh[10] = y[30]*1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  nllh[11] = y[33]*1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  nllh[12] = y[36]*1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  nllh[13] = y[39]*1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  nllh[14] = y[42]*1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  nllh[15] = 1.0/(sigmay[0]*sigmay[0]*sigmay[0])*pow(my[0]-y[0],2.0)*-1.0+1.0/sigmay[0]-y[45]*1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*5.0E-1;
  nllh[16] = y[48]*1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  nllh[17] = y[51]*1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
    break;
    case 1:
  nllh[0] = amici::log((sigmay[1]*sigmay[1])*3.141592653589793*2.0)*5.0E-1+1.0/(sigmay[1]*sigmay[1])*pow(my[1]-y[1],2.0)*5.0E-1;
  nllh[1] = y[4]*1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  nllh[2] = y[7]*1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  nllh[3] = y[10]*1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  nllh[4] = y[13]*1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  nllh[5] = y[16]*1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  nllh[6] = y[19]*1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  nllh[7] = y[22]*1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  nllh[8] = y[25]*1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  nllh[9] = y[28]*1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  nllh[10] = y[31]*1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  nllh[11] = y[34]*1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  nllh[12] = y[37]*1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  nllh[13] = y[40]*1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  nllh[14] = y[43]*1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  nllh[15] = y[46]*1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  nllh[16] = 1.0/(sigmay[1]*sigmay[1]*sigmay[1])*pow(my[1]-y[1],2.0)*-1.0+1.0/sigmay[1]-y[49]*1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*5.0E-1;
  nllh[17] = y[52]*1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
    break;
    case 2:
  nllh[0] = amici::log((sigmay[2]*sigmay[2])*3.141592653589793*2.0)*5.0E-1+1.0/(sigmay[2]*sigmay[2])*pow(my[2]-y[2],2.0)*5.0E-1;
  nllh[1] = y[5]*1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  nllh[2] = y[8]*1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  nllh[3] = y[11]*1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  nllh[4] = y[14]*1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  nllh[5] = y[17]*1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  nllh[6] = y[20]*1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  nllh[7] = y[23]*1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  nllh[8] = y[26]*1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  nllh[9] = y[29]*1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  nllh[10] = y[32]*1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  nllh[11] = y[35]*1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  nllh[12] = y[38]*1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  nllh[13] = y[41]*1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  nllh[14] = y[44]*1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  nllh[15] = y[47]*1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  nllh[16] = y[50]*1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  nllh[17] = 1.0/(sigmay[2]*sigmay[2]*sigmay[2])*pow(my[2]-y[2],2.0)*-1.0+1.0/sigmay[2]-y[53]*1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*5.0E-1;
    break;
}
}

} // namespace model_model_jakstat_adjoint_o2

} // namespace amici

