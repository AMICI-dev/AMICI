
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_jakstat_adjoint_o2{

void dJydsigma_model_jakstat_adjoint_o2(double *dJydsigma, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) {
switch(iy){
    case 0:
  dJydsigma[0+0*18] = 1.0/(sigmay[0]*sigmay[0]*sigmay[0])*pow(my[0]-y[0],2.0)*-1.0+1.0/sigmay[0];
  dJydsigma[1+0*18] = y[3]*1.0/(sigmay[0]*sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*1.0;
  dJydsigma[2+0*18] = y[6]*1.0/(sigmay[0]*sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*1.0;
  dJydsigma[3+0*18] = y[9]*1.0/(sigmay[0]*sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*1.0;
  dJydsigma[4+0*18] = y[12]*1.0/(sigmay[0]*sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*1.0;
  dJydsigma[5+0*18] = y[15]*1.0/(sigmay[0]*sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*1.0;
  dJydsigma[6+0*18] = y[18]*1.0/(sigmay[0]*sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*1.0;
  dJydsigma[7+0*18] = y[21]*1.0/(sigmay[0]*sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*1.0;
  dJydsigma[8+0*18] = y[24]*1.0/(sigmay[0]*sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*1.0;
  dJydsigma[9+0*18] = y[27]*1.0/(sigmay[0]*sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*1.0;
  dJydsigma[10+0*18] = y[30]*1.0/(sigmay[0]*sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*1.0;
  dJydsigma[11+0*18] = y[33]*1.0/(sigmay[0]*sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*1.0;
  dJydsigma[12+0*18] = y[36]*1.0/(sigmay[0]*sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*1.0;
  dJydsigma[13+0*18] = y[39]*1.0/(sigmay[0]*sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*1.0;
  dJydsigma[14+0*18] = y[42]*1.0/(sigmay[0]*sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*1.0;
  dJydsigma[15+0*18] = 1.0/(sigmay[0]*sigmay[0]*sigmay[0]*sigmay[0])*pow(my[0]-y[0],2.0)*3.0-1.0/(sigmay[0]*sigmay[0])*1.0+y[45]*1.0/(sigmay[0]*sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*1.0;
  dJydsigma[16+0*18] = y[48]*1.0/(sigmay[0]*sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*1.0;
  dJydsigma[17+0*18] = y[51]*1.0/(sigmay[0]*sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*1.0;
    break;
    case 1:
  dJydsigma[0+1*18] = 1.0/(sigmay[1]*sigmay[1]*sigmay[1])*pow(my[1]-y[1],2.0)*-1.0+1.0/sigmay[1];
  dJydsigma[1+1*18] = y[4]*1.0/(sigmay[1]*sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*1.0;
  dJydsigma[2+1*18] = y[7]*1.0/(sigmay[1]*sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*1.0;
  dJydsigma[3+1*18] = y[10]*1.0/(sigmay[1]*sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*1.0;
  dJydsigma[4+1*18] = y[13]*1.0/(sigmay[1]*sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*1.0;
  dJydsigma[5+1*18] = y[16]*1.0/(sigmay[1]*sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*1.0;
  dJydsigma[6+1*18] = y[19]*1.0/(sigmay[1]*sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*1.0;
  dJydsigma[7+1*18] = y[22]*1.0/(sigmay[1]*sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*1.0;
  dJydsigma[8+1*18] = y[25]*1.0/(sigmay[1]*sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*1.0;
  dJydsigma[9+1*18] = y[28]*1.0/(sigmay[1]*sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*1.0;
  dJydsigma[10+1*18] = y[31]*1.0/(sigmay[1]*sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*1.0;
  dJydsigma[11+1*18] = y[34]*1.0/(sigmay[1]*sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*1.0;
  dJydsigma[12+1*18] = y[37]*1.0/(sigmay[1]*sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*1.0;
  dJydsigma[13+1*18] = y[40]*1.0/(sigmay[1]*sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*1.0;
  dJydsigma[14+1*18] = y[43]*1.0/(sigmay[1]*sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*1.0;
  dJydsigma[15+1*18] = y[46]*1.0/(sigmay[1]*sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*1.0;
  dJydsigma[16+1*18] = 1.0/(sigmay[1]*sigmay[1]*sigmay[1]*sigmay[1])*pow(my[1]-y[1],2.0)*3.0-1.0/(sigmay[1]*sigmay[1])*1.0+y[49]*1.0/(sigmay[1]*sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*1.0;
  dJydsigma[17+1*18] = y[52]*1.0/(sigmay[1]*sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*1.0;
    break;
    case 2:
  dJydsigma[0+2*18] = 1.0/(sigmay[2]*sigmay[2]*sigmay[2])*pow(my[2]-y[2],2.0)*-1.0+1.0/sigmay[2];
  dJydsigma[1+2*18] = y[5]*1.0/(sigmay[2]*sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*1.0;
  dJydsigma[2+2*18] = y[8]*1.0/(sigmay[2]*sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*1.0;
  dJydsigma[3+2*18] = y[11]*1.0/(sigmay[2]*sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*1.0;
  dJydsigma[4+2*18] = y[14]*1.0/(sigmay[2]*sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*1.0;
  dJydsigma[5+2*18] = y[17]*1.0/(sigmay[2]*sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*1.0;
  dJydsigma[6+2*18] = y[20]*1.0/(sigmay[2]*sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*1.0;
  dJydsigma[7+2*18] = y[23]*1.0/(sigmay[2]*sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*1.0;
  dJydsigma[8+2*18] = y[26]*1.0/(sigmay[2]*sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*1.0;
  dJydsigma[9+2*18] = y[29]*1.0/(sigmay[2]*sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*1.0;
  dJydsigma[10+2*18] = y[32]*1.0/(sigmay[2]*sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*1.0;
  dJydsigma[11+2*18] = y[35]*1.0/(sigmay[2]*sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*1.0;
  dJydsigma[12+2*18] = y[38]*1.0/(sigmay[2]*sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*1.0;
  dJydsigma[13+2*18] = y[41]*1.0/(sigmay[2]*sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*1.0;
  dJydsigma[14+2*18] = y[44]*1.0/(sigmay[2]*sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*1.0;
  dJydsigma[15+2*18] = y[47]*1.0/(sigmay[2]*sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*1.0;
  dJydsigma[16+2*18] = y[50]*1.0/(sigmay[2]*sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*1.0;
  dJydsigma[17+2*18] = 1.0/(sigmay[2]*sigmay[2]*sigmay[2]*sigmay[2])*pow(my[2]-y[2],2.0)*3.0-1.0/(sigmay[2]*sigmay[2])*1.0+y[53]*1.0/(sigmay[2]*sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*1.0;
    break;
}
}

} // namespace model_model_jakstat_adjoint_o2

} // namespace amici

