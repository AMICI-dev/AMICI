#ifndef amici_symbolic_functions_h
#define amici_symbolic_functions_h

namespace amici {

double amilog(double x);
double dirac(double x);
double heaviside(double x);

double am_min(double a, double b, double c);
double Dam_min(int id, double a, double b, double c);
double am_max(double a, double b, double c);
double Dam_max(int id, double a, double b, double c);

int amiIsNaN(double what);
int amiIsInf(double what);
double amiGetNaN();

/* sign */
double sign(double x);

/* splines */

double am_spline(double t, int num, ...);
double am_spline_pos(double t, int num, ...);
double am_Dspline(int id, double t, int num, ...);
double am_Dspline_pos(int id, double t, int num, ...);
double am_DDspline(int id1, int id2, double t, int num, ...);
double am_DDspline_pos(int id1, int id2, double t, int num, ...);
} // namespace amici

#endif /* amici_symbolic_functions_h */
