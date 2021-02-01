/**
 * @file   symbolic_functions.cpp
 * @brief  definition of symbolic functions
 *
 * This file contains definitions of various symbolic functions which
 */

#include "amici/symbolic_functions.h"
#include "amici/spline.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdarg>
#include <cstdlib>
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

/**
 * c++ interface to the isNaN function
 *
 * @param what argument
 * @return isnan(what)
 *
 */
int isNaN(double what) { return std::isnan(what); }

/**
 * c++ interface to the isinf function
 *
 * @param what argument
 * @return isnan(what)
 *
 */
int isInf(double what) { return std::isinf(what); }

/**
 * function returning nan
 *
 * @return NaN
 *
 */
double getNaN() { return NAN; }

/**
 * C implementation of log function, this prevents returning NaN values for
 * negative values
 *
 * @param x argument
 * @return if(x>0) then log(x) else -Inf
 *
 */
double log(double x) {
    if (x <= 0) {
        return -std::log(DBL_MAX);
    }
    return std::log(x);
}

/**
 * C implementation of matlab function dirac
 *
 * @param x argument
 * @return if(x==0) then INF else 0
 *
 */
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
 * @return if(x>0) then 1 else 0
 *
 */
double heaviside(double x) {
    if (x < 0.0) {
        return 0.0;
    }
    return 1.0;
}

/**
 *  c implementation of matlab function sign
 *
 * @param x argument
 * @return 0
 *
 */
double sign(double x) {
    if (x > 0.0)
        return 1.0;

    if (x < 0.0)
        return -1.0;

    return 0.0;
}

/**
 * c implementation of matlab function max
 *
 * @param a value1
 * @param b value2
 * @return if(a > b) then a else b
 *
 */
double max(double a, double b, double  /*c*/) {
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

/**
 * c implementation of matlab function min
 *
 * @param a value1
 * @param b value2
 * @param c bogus parameter to ensure correct parsing as a function
 * @return if(a < b) then a else b
 *
 */
double min(double a, double b, double c) {
    return (-max(-a,-b,c));
}

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
double Dmax(int id, double a, double b, double  /*c*/) {
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
double Dmin(int id, double a, double b, double c) {
    return Dmax(id,-a,-b,c);
}

/**
 * specialized pow functions that assumes positivity of the first argument
 *
 * @param base base
 * @param exponent exponent
 * @return pow(max(base,0.0),exponent)
 *
 */
double pos_pow(double base, double exponent) {
    // we do NOT want to propagate NaN values here, if base is nan, so should the output be
    return pow(std::max(base, 0.0),exponent);
}

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
double spline(double t, int num, ...) {

    va_list valist;

    double uout;
    double ss;
    double dudt;

    auto *ts = (double *)alloca(num * sizeof(double));
    auto *us = (double *)alloca(num * sizeof(double));

    auto *b = (double *)alloca(num * sizeof(double));
    auto *c = (double *)alloca(num * sizeof(double));
    auto *d = (double *)alloca(num * sizeof(double));

    /* Variable list type macro */
    /* initialize valist for num number of arguments */
    va_start(valist, num);

    for (int i = 0; i < 2 * num; i += 2) {
        int j = i / 2;
        ts[j] = va_arg(valist, double);
        us[j] = va_arg(valist, double);
    }
    ss = va_arg(valist, double);
    dudt = va_arg(valist, double);

    /* clean memory reserved for valist */
    va_end(valist);

    spline(num, static_cast<int>(ss), 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(num, t, ts, us, b, c, d);

    return uout;
}

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
double spline_pos(double t, int num, ...) {

    va_list valist;

    double uout;
    double ss;
    double dudt;

    auto *ts = (double *)alloca(num * sizeof(double));
    auto *us = (double *)alloca(num * sizeof(double));
    auto *uslog = (double *)alloca(num * sizeof(double));

    auto *b = (double *)alloca(num * sizeof(double));
    auto *c = (double *)alloca(num * sizeof(double));
    auto *d = (double *)alloca(num * sizeof(double));

    /* initialize valist for num number of arguments */
    va_start(valist, num);

    for (int i = 0; i < 2 * num; i += 2) {
        int j = i / 2;
        ts[j] = va_arg(valist, double);
        us[j] = va_arg(valist, double);
        uslog[j] = log(us[j]);
    }
    ss = va_arg(valist, double);
    dudt = va_arg(valist, double);

    /* clean memory reserved for valist */
    va_end(valist);

    spline(num, static_cast<int>(ss), 0, dudt, 0.0, ts, uslog, b, c, d);
    uout = seval(num, t, ts, uslog, b, c, d);

    return exp(uout);
}

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
double Dspline(int id, double t, int num, ...) {

    va_list valist;

    double uout;
    double ss;
    double dudt;

    double *ts = (double *)alloca(num * sizeof(double));
    double *us = (double *)alloca(num * sizeof(double));
    double *ps = (double *)alloca(num * sizeof(double));

    double *b = (double *)alloca(num * sizeof(double));
    double *c = (double *)alloca(num * sizeof(double));
    double *d = (double *)alloca(num * sizeof(double));

    int did = id / 2 - 2;

    /* initialize valist for num number of arguments */
    va_start(valist, num);

    for (int i = 0; i < 2 * num; i += 2) {
        int j = i / 2;
        ts[j] = va_arg(valist, double);
        ps[j] = va_arg(valist, double);
        us[j] = 0.0;
    }
    us[did] = 1.0;
    ss = va_arg(valist, double);
    dudt = va_arg(valist, double);

    /* clean memory reserved for valist */
    va_end(valist);

    spline(num, static_cast<int>(ss), 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(num, t, ts, us, b, c, d);

    return uout;
}

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
double Dspline_pos(int id, double t, int num, ...) {

    va_list valist;

    auto *ts = (double *)alloca(num * sizeof(double));
    auto *us = (double *)alloca(num * sizeof(double));
    auto *sus = (double *)alloca(num * sizeof(double));
    auto *uslog = (double *)alloca(num * sizeof(double));

    auto *b = (double *)alloca(num * sizeof(double));
    auto *c = (double *)alloca(num * sizeof(double));
    auto *d = (double *)alloca(num * sizeof(double));

    double uout;
    double ss;
    double dudt;
    double uspline_pos;
    double suspline;

    int did = id / 2 - 2;

    /* initialize valist for num number of arguments */
    va_start(valist, num);

    for (int i = 0; i < 2 * num; i += 2) {
        int j = i / 2;
        ts[j] = va_arg(valist, double);
        us[j] = va_arg(valist, double);
        uslog[j] = log(us[j]);
        sus[j] = 0.0;
    }
    ss = va_arg(valist, double);
    dudt = va_arg(valist, double);
    sus[did] = 1.0;

    /* clean memory reserved for valist */
    va_end(valist);

    spline(num, static_cast<int>(ss), 0, dudt, 0.0, ts, uslog, b, c, d);
    uspline_pos = exp(seval(num, t, ts, uslog, b, c, d));

    spline(num, static_cast<int>(ss), 0, dudt, 0.0, ts, sus, b, c, d);
    suspline = seval(num, t, ts, sus, b, c, d);
    uout = suspline * uspline_pos / us[did];

    return uout;
}

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
double DDspline(int  /*id1*/, int  /*id2*/, double  /*t*/, int  /*num*/, ...) { return 0.0; }

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
double DDspline_pos(int id1, int id2, double t, int num, ...) {

    va_list valist;

    auto *ts = (double *)alloca(num * sizeof(double));
    auto *us = (double *)alloca(num * sizeof(double));
    auto *sus1 = (double *)alloca(num * sizeof(double));
    auto *sus2 = (double *)alloca(num * sizeof(double));
    auto *uslog = (double *)alloca(num * sizeof(double));

    auto *b = (double *)alloca(num * sizeof(double));
    auto *c = (double *)alloca(num * sizeof(double));
    auto *d = (double *)alloca(num * sizeof(double));

    double uout;
    double ss;
    double dudt;
    double uspline_pos;
    double su1spline;
    double su2spline;

    int did1 = id1 / 2 - 2;
    int did2 = id2 / 2 - 2;

    /* initialize valist for num number of arguments */
    va_start(valist, num);

    for (int i = 0; i < 2 * num; i += 2) {
        int j = i / 2;
        ts[j] = va_arg(valist, double);
        us[j] = va_arg(valist, double);
        uslog[j] = log(us[j]);
        sus1[j] = 0.0;
        sus2[j] = 0.0;
    }
    ss = va_arg(valist, double);
    dudt = va_arg(valist, double);
    sus1[did1] = 1.0;
    sus2[did2] = 1.0;

    /* clean memory reserved for valist */
    va_end(valist);

    spline(num, static_cast<int>(ss), 0, dudt, 0.0, ts, uslog, b, c, d);
    uspline_pos = exp(seval(num, t, ts, uslog, b, c, d));

    spline(num, static_cast<int>(ss), 0, dudt, 0.0, ts, sus1, b, c, d);
    su1spline = seval(num, t, ts, sus1, b, c, d);

    spline(num, static_cast<int>(ss), 0, dudt, 0.0, ts, sus2, b, c, d);
    su2spline = seval(num, t, ts, sus2, b, c, d);

    if (id1 == id2) {
        uout = (su1spline * su2spline - su1spline) * uspline_pos;
    } else {
        uout = su1spline * su2spline * uspline_pos;
    }
    uout = uout / us[did1] / us[did2];

    return uout;
}

} // namespace amici
