/**
 * @file   symbolic_functions.cpp
 * @brief  definition of symbolic functions
 *
 * This file contains definitions of various symbolic functions which
 */

#include "amici/symbolic_functions.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace amici {

double get_nan() { return std::numeric_limits<double>::quiet_NaN(); }

double log(double x) {
    if (x <= 0) {
        return -std::log(std::numeric_limits<double>::max());
    }
    return std::log(x);
}

double dirac(double x) {
    if (x == 0.0) {
        return std::numeric_limits<double>::max();
    }
    return 0.0;
}

double heaviside(double x, double x0) {
    if (x < 0.0)
        return 0.0;
    if (x == 0.0)
        return x0;
    return 1.0;
}

double sign(double x) {
    if (x > 0.0)
        return 1.0;

    if (x < 0.0)
        return -1.0;

    return 0.0;
}

double pos_pow(double base, double exponent) {
    // we do NOT want to propagate NaN values here, if base is NaN, so should
    // the output be
    return pow(std::max(base, 0.0), exponent);
}

} // namespace amici
