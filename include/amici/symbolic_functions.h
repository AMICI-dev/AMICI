#ifndef amici_symbolic_functions_h
#define amici_symbolic_functions_h

namespace amici {

/**
 * C implementation of log function, this prevents returning NaN values for
 * negative values
 *
 * @param x argument
 * @return if(x>0) then log(x) else -Inf
 *
 */
double log(double x);

/**
 * C implementation of matlab function dirac
 *
 * @param x argument
 * @return if(x==0) then INF else 0
 *
 */
double dirac(double x);
double heaviside(double x, double x0);

/**
 * c implementation of matlab function min
 *
 * @param a value1
 * @param b value2
 * @param c bogus parameter to ensure correct parsing as a function
 * @return if(a < b) then a else b
 *
 */
double min(double a, double b, double c);

/**
 * parameter derivative of c implementation of matlab function max
 *
 * @param id argument index for differentiation
 * @param a value1
 * @param b value2
 * @param c bogus parameter to ensure correct parsing as a function
 * @return id == 1:  if(a > b) then 1 else 0
 * @return id == 2:  if(a > b) then 0 else 1
 *
 */
double Dmin(int id, double a, double b, double c);

/**
 * c implementation of matlab function max
 *
 * @param a value1
 * @param b value2
 * @return if(a > b) then a else b
 *
 */
double max(double a, double b, double c);

/**
 * parameter derivative of c implementation of matlab function max
 *
 * @param id argument index for differentiation
 * @param a value1
 * @param b value2
 * @return id == 1:  if(a > b) then 1 else 0
 * @return id == 2:  if(a > b) then 0 else 1
 *
 */
double Dmax(int id, double a, double b, double c);

/**
 * specialized pow functions that assumes positivity of the first argument
 *
 * @param base base
 * @param exponent exponent
 * @return pow(max(base,0.0),exponen)
 *
 */
double pos_pow(double base, double exponent);

/**
 * c++ interface to the isNaN function
 *
 * @param what argument
 * @return isnan(what)
 *
 */
int isNaN(double what);

/**
 * c++ interface to the isinf function
 *
 * @param what argument
 * @return isnan(what)
 *
 */
int isInf(double what);

/**
 * function returning nan
 *
 * @return NaN
 *
 */
double getNaN();

/**
 *  c implementation of matlab function sign
 *
 * @param x argument
 * @return 0
 *
 */
double sign(double x);

} // namespace amici

#endif /* amici_symbolic_functions_h */
