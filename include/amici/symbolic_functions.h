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

/* legacy spline implementation in C (MATLAB only) */

/**
 * @brief Spline function
 *
 * Takes variable argument pairs (ti,pi) with `ti`: location of node i and
 * `pi`: spline value at node i. the last two arguments are always `ss`: flag
 * indicating whether slope at first node should be user defined
 * and `dudt` user defined slope at first node. All arguments must be of type
 * double.
 *
 * @param t point at which the spline should be evaluated
 * @param num number of spline nodes
 *
 * @return spline(t)
 *
 */
double spline(double t, int num, ...);

/**
 * @brief Exponentiated spline function
 *
 * Takes variable argument pairs (ti,pi) with `ti`: location of node i and
 * `pi`: spline value at node i. the last two arguments are always `ss`: flag
 * indicating whether slope at first node should be user defined
 * and `dudt` user defined slope at first node. All arguments must be of type
 * double.
 *
 * @param t point at which the spline should be evaluated
 * @param num number of spline nodes
 *
 * @return spline(t)
 *
 */
double spline_pos(double t, int num, ...);

/**
 * @brief Derivation of a spline function
 *
 * Takes variable argument pairs (ti,pi) with `ti`: location of node i and
 * `pi`: spline value at node i. the last two arguments are always `ss`: flag
 * indicating whether slope at first node should be user defined
 * and `dudt` user defined slope at first node. All arguments but id must be of
 * type double.
 *
 * @param id index of node to which the derivative of the corresponding spline
 * coefficient should be computed
 * @param t point at which the spline should be evaluated
 * @param num number of spline nodes
 *
 * @return dsplinedp(t)
 *
 */
double Dspline(int id, double t, int num, ...);

/**
 * @brief Derivation of an exponentiated spline function
 *
 * Takes variable argument pairs (ti,pi) with `ti`: location of node i and
 * `pi`: spline value at node i. the last two arguments are always `ss`: flag
 * indicating whether slope at first node should be user defined
 * and `dudt` user defined slope at first node. All arguments but id must be of
 * type double.
 *
 * @param id index of node to which the derivative of the corresponding spline
 * coefficient should be computed
 * @param t point at which the spline should be evaluated
 * @param num number of spline nodes
 *
 * @return dsplinedp(t)
 *
 */
double Dspline_pos(int id, double t, int num, ...);

/**
 * @brief Second derivation of a spline function
 *
 * Takes variable argument pairs (ti,pi) with `ti`: location of node i and
 * `pi`: spline value at node i. the last two arguments are always `ss`: flag
 * indicating whether slope at first node should be user defined
 * and `dudt` user defined slope at first node. All arguments but id1 and id2
 * must be of type double.
 *
 * @param id1 index of node to which the first derivative of the corresponding
 * spline coefficient should be computed
 * @param id2 index of node to which the second derivative of the corresponding
 * spline coefficient should be computed
 * @param t point at which the spline should be evaluated
 * @param num number of spline nodes
 *
 * @return ddspline(t)
 */
double DDspline(int id1, int id2, double t, int num, ...);

/**
 * @brief Derivation of an exponentiated spline function
 *
 * Takes variable argument pairs (ti,pi) with `ti`: location of node i and
 * `pi`: spline value at node i. the last two arguments are always `ss`: flag
 * indicating whether slope at first node should be user defined
 * and `dudt` user defined slope at first node. All arguments but id1 and id2
 * must be of type double.
 *
 * @param id1 index of node to which the first derivative of the corresponding
 * spline coefficient should be computed
 * @param id2 index of node to which the second derivative of the corresponding
 * spline coefficient should be computed
 * @param t point at which the spline should be evaluated
 * @param num number of spline nodes
 *
 * @return ddspline(t)
 *
 */
double DDspline_pos(int id1, int id2, double t, int num, ...);

} // namespace amici

#endif /* amici_symbolic_functions_h */
