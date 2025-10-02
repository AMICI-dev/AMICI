/**
 * @file   symbolic_functions.cpp
 * @brief  definition of symbolic functions
 *
 * This file contains definitions of various symbolic functions which
 */

#include "amici/symbolic_functions.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdarg>
#if _MSC_VER && !__INTEL_COMPILER
#include <malloc.h>
#define alloca _alloca
#elif defined(WIN32) || defined(__WIN32) || defined(__WIN32__)
// For alloca().
#include <malloc.h>
#else
#include <alloca.h>
#endif

namespace amici {

int isNaN(double what) { return std::isnan(what); }

int isInf(double what) { return std::isinf(what); }

double getNaN() { return NAN; }

double log(double x) {
    if (x <= 0) {
        return -std::log(DBL_MAX);
    }
    return std::log(x);
}

double dirac(double x) {
    if (x == 0.0) {
        return DBL_MAX;
    }
    return 0.0;
}

/**
 * c implementation of matlab function heaviside
 *
 * @param x argument
 * @param x0 value at x==0
 * @return if(x>0) then 1 else if (x==0) then x0 else 0
 *
 */
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

double max(double a, double b, double /*c*/) {
    int anan = isNaN(a), bnan = isNaN(b);
    if (anan || bnan) {
        if (anan && !bnan)
            return b;
        if (!anan && bnan)
            return a;
        return a;
    }
    return (std::max(a, b));
}

double min(double a, double b, double c) { return (-max(-a, -b, c)); }

double Dmax(int id, double a, double b, double /*c*/) {
    if (id == 1.0) {
        if (a > b)
            return 1.0;
        return 0.0;
    }
    if (a > b) {
        return 0.0;
    }
    return 1.0;
}

double Dmin(int id, double a, double b, double c) {
    return Dmax(id, -a, -b, c);
}

double pos_pow(double base, double exponent) {
    // we do NOT want to propagate NaN values here, if base is nan, so should
    // the output be
    return pow(std::max(base, 0.0), exponent);
}

} // namespace amici
