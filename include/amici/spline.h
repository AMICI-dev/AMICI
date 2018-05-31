#ifndef amici_spline_h
#define amici_spline_h
#include <math.h>

namespace amici {

int spline(int n, int end1, int end2, double slope1, double slope2, double x[],
           double y[], double b[], double c[], double d[]);

double seval(int n, double u, double x[], double y[], double b[], double c[],
             double d[]);

double sinteg(int n, double u, double x[], double y[], double b[], double c[],
              double d[]);

} // namespace amici

#endif /* spline_h */
