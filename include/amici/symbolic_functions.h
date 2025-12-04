#ifndef AMICI_SYMBOLIC_FUNCTIONS_H
#define AMICI_SYMBOLIC_FUNCTIONS_H

namespace amici {

/**
 * @brief A log function, that prevents returning NaN values for
 * negative values.
 *
 * @param x argument
 * @return if(x>0) then log(x) else -Inf
 */
double log(double x);

/**
 * @brief The Dirac delta function.
 *
 * @param x argument
 * @return if(x==0) then INF else 0
 *
 */
double dirac(double x);

/**
 * @brief The Heaviside function.
 *
 * @param x argument
 * @param x0 value at x==0
 * @return if(x>0) then 1 else if (x==0) then x0 else 0
 *
 */
double heaviside(double x, double x0);

/**
 * @brief Specialized pow functions that assumes positivity of the first
 * argument.
 *
 * @param base base
 * @param exponent exponent
 * @return pow(max(base,0.0), exponent)
 *
 */
double pos_pow(double base, double exponent);

/**
 * @brief Get NaN value.
 *
 * @return NaN
 */
double get_nan();

/**
 * @brief The sign function.
 *
 * @param x argument
 * @return if(x>0) then 1 else if (x<0) then -1 else 0
 *
 */
double sign(double x);

} // namespace amici

#endif // AMICI_SYMBOLIC_FUNCTIONS_H
