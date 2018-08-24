#ifndef amici_symbolic_functions_h
#define amici_symbolic_functions_h

namespace amici {

double log(double x);
double dirac(double x);
double heaviside(double x);

double min(double a, double b, double c);
double Dmin(int id, double a, double b, double c);
double max(double a, double b, double c);
double Dmax(int id, double a, double b, double c);

double pos_pow(double base, double exponent);
    
int isNaN(double what);
int isInf(double what);
double getNaN();

/* sign */
double sign(double x);

/* splines */

double spline(double t, int num, ...);
double spline_pos(double t, int num, ...);
double Dspline(int id, double t, int num, ...);
double Dspline_pos(int id, double t, int num, ...);
double DDspline(int id1, int id2, double t, int num, ...);
double DDspline_pos(int id1, int id2, double t, int num, ...);
} // namespace amici

#endif /* amici_symbolic_functions_h */
