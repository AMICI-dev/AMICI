
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_jakstat_adjoint_o2{

void dJydy_model_jakstat_adjoint_o2(double *dJydy, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) {
switch(iy){
    case 0:
  dJydy[0+0*18] = 1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  dJydy[1+0*18] = y[3]*1.0/(sigmay[0]*sigmay[0])*1.0;
  dJydy[1+3*18] = 1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  dJydy[2+0*18] = y[6]*1.0/(sigmay[0]*sigmay[0])*1.0;
  dJydy[2+6*18] = 1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  dJydy[3+0*18] = y[9]*1.0/(sigmay[0]*sigmay[0])*1.0;
  dJydy[3+9*18] = 1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  dJydy[4+0*18] = y[12]*1.0/(sigmay[0]*sigmay[0])*1.0;
  dJydy[4+12*18] = 1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  dJydy[5+0*18] = y[15]*1.0/(sigmay[0]*sigmay[0])*1.0;
  dJydy[5+15*18] = 1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  dJydy[6+0*18] = y[18]*1.0/(sigmay[0]*sigmay[0])*1.0;
  dJydy[6+18*18] = 1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  dJydy[7+0*18] = y[21]*1.0/(sigmay[0]*sigmay[0])*1.0;
  dJydy[7+21*18] = 1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  dJydy[8+0*18] = y[24]*1.0/(sigmay[0]*sigmay[0])*1.0;
  dJydy[8+24*18] = 1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  dJydy[9+0*18] = y[27]*1.0/(sigmay[0]*sigmay[0])*1.0;
  dJydy[9+27*18] = 1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  dJydy[10+0*18] = y[30]*1.0/(sigmay[0]*sigmay[0])*1.0;
  dJydy[10+30*18] = 1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  dJydy[11+0*18] = y[33]*1.0/(sigmay[0]*sigmay[0])*1.0;
  dJydy[11+33*18] = 1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  dJydy[12+0*18] = y[36]*1.0/(sigmay[0]*sigmay[0])*1.0;
  dJydy[12+36*18] = 1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  dJydy[13+0*18] = y[39]*1.0/(sigmay[0]*sigmay[0])*1.0;
  dJydy[13+39*18] = 1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  dJydy[14+0*18] = y[42]*1.0/(sigmay[0]*sigmay[0])*1.0;
  dJydy[14+42*18] = 1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  dJydy[15+0*18] = 1.0/(sigmay[0]*sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*1.0+y[45]*1.0/(sigmay[0]*sigmay[0])*1.0;
  dJydy[15+45*18] = 1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  dJydy[16+0*18] = y[48]*1.0/(sigmay[0]*sigmay[0])*1.0;
  dJydy[16+48*18] = 1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  dJydy[17+0*18] = y[51]*1.0/(sigmay[0]*sigmay[0])*1.0;
  dJydy[17+51*18] = 1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
    break;
    case 1:
  dJydy[0+1*18] = 1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  dJydy[1+1*18] = y[4]*1.0/(sigmay[1]*sigmay[1])*1.0;
  dJydy[1+4*18] = 1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  dJydy[2+1*18] = y[7]*1.0/(sigmay[1]*sigmay[1])*1.0;
  dJydy[2+7*18] = 1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  dJydy[3+1*18] = y[10]*1.0/(sigmay[1]*sigmay[1])*1.0;
  dJydy[3+10*18] = 1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  dJydy[4+1*18] = y[13]*1.0/(sigmay[1]*sigmay[1])*1.0;
  dJydy[4+13*18] = 1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  dJydy[5+1*18] = y[16]*1.0/(sigmay[1]*sigmay[1])*1.0;
  dJydy[5+16*18] = 1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  dJydy[6+1*18] = y[19]*1.0/(sigmay[1]*sigmay[1])*1.0;
  dJydy[6+19*18] = 1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  dJydy[7+1*18] = y[22]*1.0/(sigmay[1]*sigmay[1])*1.0;
  dJydy[7+22*18] = 1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  dJydy[8+1*18] = y[25]*1.0/(sigmay[1]*sigmay[1])*1.0;
  dJydy[8+25*18] = 1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  dJydy[9+1*18] = y[28]*1.0/(sigmay[1]*sigmay[1])*1.0;
  dJydy[9+28*18] = 1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  dJydy[10+1*18] = y[31]*1.0/(sigmay[1]*sigmay[1])*1.0;
  dJydy[10+31*18] = 1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  dJydy[11+1*18] = y[34]*1.0/(sigmay[1]*sigmay[1])*1.0;
  dJydy[11+34*18] = 1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  dJydy[12+1*18] = y[37]*1.0/(sigmay[1]*sigmay[1])*1.0;
  dJydy[12+37*18] = 1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  dJydy[13+1*18] = y[40]*1.0/(sigmay[1]*sigmay[1])*1.0;
  dJydy[13+40*18] = 1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  dJydy[14+1*18] = y[43]*1.0/(sigmay[1]*sigmay[1])*1.0;
  dJydy[14+43*18] = 1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  dJydy[15+1*18] = y[46]*1.0/(sigmay[1]*sigmay[1])*1.0;
  dJydy[15+46*18] = 1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  dJydy[16+1*18] = 1.0/(sigmay[1]*sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*1.0+y[49]*1.0/(sigmay[1]*sigmay[1])*1.0;
  dJydy[16+49*18] = 1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
  dJydy[17+1*18] = y[52]*1.0/(sigmay[1]*sigmay[1])*1.0;
  dJydy[17+52*18] = 1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
    break;
    case 2:
  dJydy[0+2*18] = 1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  dJydy[1+2*18] = y[5]*1.0/(sigmay[2]*sigmay[2])*1.0;
  dJydy[1+5*18] = 1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  dJydy[2+2*18] = y[8]*1.0/(sigmay[2]*sigmay[2])*1.0;
  dJydy[2+8*18] = 1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  dJydy[3+2*18] = y[11]*1.0/(sigmay[2]*sigmay[2])*1.0;
  dJydy[3+11*18] = 1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  dJydy[4+2*18] = y[14]*1.0/(sigmay[2]*sigmay[2])*1.0;
  dJydy[4+14*18] = 1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  dJydy[5+2*18] = y[17]*1.0/(sigmay[2]*sigmay[2])*1.0;
  dJydy[5+17*18] = 1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  dJydy[6+2*18] = y[20]*1.0/(sigmay[2]*sigmay[2])*1.0;
  dJydy[6+20*18] = 1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  dJydy[7+2*18] = y[23]*1.0/(sigmay[2]*sigmay[2])*1.0;
  dJydy[7+23*18] = 1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  dJydy[8+2*18] = y[26]*1.0/(sigmay[2]*sigmay[2])*1.0;
  dJydy[8+26*18] = 1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  dJydy[9+2*18] = y[29]*1.0/(sigmay[2]*sigmay[2])*1.0;
  dJydy[9+29*18] = 1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  dJydy[10+2*18] = y[32]*1.0/(sigmay[2]*sigmay[2])*1.0;
  dJydy[10+32*18] = 1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  dJydy[11+2*18] = y[35]*1.0/(sigmay[2]*sigmay[2])*1.0;
  dJydy[11+35*18] = 1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  dJydy[12+2*18] = y[38]*1.0/(sigmay[2]*sigmay[2])*1.0;
  dJydy[12+38*18] = 1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  dJydy[13+2*18] = y[41]*1.0/(sigmay[2]*sigmay[2])*1.0;
  dJydy[13+41*18] = 1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  dJydy[14+2*18] = y[44]*1.0/(sigmay[2]*sigmay[2])*1.0;
  dJydy[14+44*18] = 1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  dJydy[15+2*18] = y[47]*1.0/(sigmay[2]*sigmay[2])*1.0;
  dJydy[15+47*18] = 1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  dJydy[16+2*18] = y[50]*1.0/(sigmay[2]*sigmay[2])*1.0;
  dJydy[16+50*18] = 1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
  dJydy[17+2*18] = 1.0/(sigmay[2]*sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*1.0+y[53]*1.0/(sigmay[2]*sigmay[2])*1.0;
  dJydy[17+53*18] = 1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
    break;
}
}

} // namespace model_model_jakstat_adjoint_o2

} // namespace amici

