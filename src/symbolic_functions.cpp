/**
 * @file   symbolic_functions.cpp
 * @brief  definition of symbolic functions
 *
 * This file contains definitions of various symbolic functions which
 */


#include <cstdarg>
#include <algorithm>
#include <cmath>
#ifndef AMICI_WITHOUT_MATLAB
    #include <mex.h>
#endif
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cfloat>
#include <alloca.h>
#include <include/symbolic_functions.h>
#include <include/spline.h>


#undef ts

/*! bool return value true */
#define TRUE 1
/*! bool return value false */
#define FALSE 0


int amiIsNaN(double what) {
#if defined mex_typedefs_h || defined mex_h
    return mxIsNaN(what);
#else
    return std::isnan(what);
#endif
}

int amiIsInf(double what) {
#if defined mex_typedefs_h || defined mex_h
    return mxIsInf(what);
#else
    return std::isinf(what);
#endif
}

double amiGetNaN() {
#if defined mex_typedefs_h || defined mex_h
    return mxGetNaN();
#else
    return INFINITY;
#endif
}

void fillArray(double *destination, int count, double value) {
    int i;
    for(i = 0; i < count; ++i)
        destination[i] = value;
}

double sum(double const *array, int numElements) {
    double sum = 0;
    int i;
    for(i = 0; i < numElements; ++i) {
        sum += array[i];
    }
    return sum;
}

void zeros(double *destination, int count) {
    memset(destination, 0, sizeof(double) * count);
}

void ones(double *destination, int count) {
    fillArray(destination, count, 1);
}

void linSpace(double *destination, double from, double to, int numValues) {
    double delta = (to - from) / (numValues - 1);
    int i;
    for(i = 0; i < numValues; ++i) {
        destination[i] = from + i * delta;
    }
}

double *linSpaceAlloc(double from, double to, int numValues) {
    double *destination = (double*)malloc(sizeof(double) * numValues);
    linSpace(destination, from, to, numValues);
    return destination;
}

void printArray(double const *array, int numElements) {
    printfArray(array, numElements, "%e\t");
}

void printfArray(double const *array, int numElements, char const *format) {
    int i;
    for(i = 0; i < numElements; ++i) {
        printf(format, array[i]);
    }
}

void errMsgIdAndTxt(
    const char * identifier, /* string with error message identifier */
    const char * err_msg,    /* string with error message printf-style format */
    ...                      /* any additional arguments */
    ) {
#ifdef AMICI_WITHOUT_MATLAB
    printf("[Error] %s: %s\n", identifier, err_msg);
#else
    mexWarnMsgIdAndTxt(identifier, err_msg);
#endif
}

void warnMsgIdAndTxt(
    const char * identifier, /* string with error message identifier */
    const char * err_msg,    /* string with error message printf-style format */
    ...                      /* any additional arguments */
    ) {
#ifdef AMICI_WITHOUT_MATLAB
    printf("[Warning] %s: %s\n", identifier, err_msg);
#else
    mexWarnMsgIdAndTxt(identifier, err_msg);
#endif
}

/**
 * c implementation of log function, this prevents returning NaN values for negative values
 *
 * @param x argument
 * @return if(x>0) then log(x) else -Inf
 *
 */
double amilog(double x) {
    if (x<=0) {
        return(-log(DBL_MAX));
    } else {
        return(log(x));
    }
}

/**
 * c implementation of matlab function dirac
 *
 * @param x argument
 * @return if(x==0) then INF else 0
 *
 */
double dirac(double x) {
    if (x == 0) {
        return(DBL_MAX);
    } else {
        return(0);
    }
}

/**
 * c implementation of matlab function heaviside
 *
 * @param x argument
 * @return if(x>0) then 1 else 0
 *
 */
double heaviside(double x) {
    if (x <= 0) {
        return(0);
    } else {
        return(1);
    }
}



/**
 *  c implementation of matlab function sign
 *
 * @param x argument
 * @return 0 @type double
 *
 */
double sign(double x) {
    if (x > 0) {
        return(1);
    } else {
        if (x < 0 ) {
            return(-1);
        } else {
            return(0);
        }
    }
}

/**
 * c implementation of matlab function min
 *
 * @param a value1 @type double
 * @param b value2 @type double
 * @param c bogus parameter to ensure correct parsing as a function @type double
 * @return if(a < b) then a else b @type double
 *
 */
double am_min(double a, double b, double c) {
    int anan = amiIsNaN(a), bnan = amiIsNaN(b);
    if(anan || bnan) {
        if(anan && !bnan) return b;
        if(!anan && bnan) return a;
        return a;
    }
    return(std::min(a,b));
}

/**
 * parameter derivative of c implementation of matlab function min
 *
 * @param id argument index for differentiation
 * @param a value1 @type double
 * @param b value2 @type double
 * @param c bogus parameter to ensure correct parsing as a function @type double
 * @return id == 1:  if(a < b) then 1 else 0 @type double
 * @return id == 2:  if(a < b) then 0 else 1 @type double
 *
 */
double Dam_min(int id,double a, double b, double c) {
    if (id == 1) {
        if (a < b) {
            return(1);
        } else {
            return(0);
        }
    } else {
        if (a < b) {
            return(0);
        } else {
            return(1);
        }
    }
}

/**
 * c implementation of matlab function max
 *
 * @param a value1 @type double
 * @param b value2 @type double
 * @param c bogus parameter to ensure correct parsing as a function @type double
 * @return if(a > b) then a else b @type double
 *
 */
double am_max(double a, double b, double c) {
    int anan = amiIsNaN(a), bnan = amiIsNaN(b);
    if(anan || bnan) {
        if(anan && !bnan) return b;
        if(!anan && bnan) return a;
        return a;
    }
    return(std::max(a,b));
}

/**
 * parameter derivative of c implementation of matlab function max
 *
 * @param id argument index for differentiation
 * @param a value1 @type double
 * @param b value2 @type double
 * @param c bogus parameter to ensure correct parsing as a function @type double
 * @return id == 1:  if(a > b) then 1 else 0 @type double
 * @return id == 2:  if(a > b) then 0 else 1 @type double
 *
 */
double Dam_max(int id,double a, double b, double c) {
    if (id == 1) {
        if (a > b) {
            return(1);
        } else {
            return(0);
        }
    } else {
        if (a > b) {
            return(0);
        } else {
            return(1);
        }
    }
}


/**
 * spline function
 *
 * @param t point at which the spline should be evaluated
 * @param ti location of node i
 * @param pi spline value at node i
 * @param ss flag indicating whether slope at first node should be user defined
 * @param dudt user defined slope at first node
 *
 * @return spline(t)
 *
 */
double am_spline(double t, int num, ...) {
    
    va_list valist;
    
    double uout;
    double ss;
    double dudt;
    
    double *ts = (double*) alloca(num*sizeof(double));
    double *us = (double*) alloca(num*sizeof(double));
    
    double *b = (double*) alloca(num*sizeof(double));
    double *c = (double*) alloca(num*sizeof(double));
    double *d = (double*) alloca(num*sizeof(double));
    
    int i;
    int j;
    
     /* Variable list type macro */
    /* initialize valist for num number of arguments */
    va_start(valist, num);
    
    for (i=0; i<2*num; i+=2) {
        j = i/2;
        ts[j] = va_arg(valist, double);
        us[j] = va_arg(valist, double);
    }
    ss = va_arg(valist, double);
    dudt = va_arg(valist, double);
    
    /* clean memory reserved for valist */
    va_end(valist);
    
    spline(num, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(num, t, ts, us, b, c, d);
    
    return(uout);
}

/**
 * exponentiated spline function
 *
 * @param t point at which the spline should be evaluated
 * @param ti location of node i
 * @param pi spline value at node i
 * @param ss flag indicating whether slope at first node should be user defined
 * @param dudt user defined slope at first node
 *
 * @return spline(t)
 *
 */
double am_spline_pos(double t, int num, ...) {
    
    va_list valist;
    
    double uout;
    double ss;
    double dudt;
    
    double *ts = (double*) alloca(num*sizeof(double));
    double *us = (double*) alloca(num*sizeof(double));
    double *uslog = (double*) alloca(num*sizeof(double));
    
    double *b = (double*) alloca(num*sizeof(double));
    double *c = (double*) alloca(num*sizeof(double));
    double *d = (double*) alloca(num*sizeof(double));
    
    int i;
    int j;
    
    /* initialize valist for num number of arguments */
    va_start(valist, num);
    
    for (i=0; i<2*num; i+=2) {
        j = i/2;
        ts[j] = va_arg(valist, double);
        us[j] = va_arg(valist, double);
        uslog[j] = log(us[j]);
    }
    ss = va_arg(valist, double);
    dudt = va_arg(valist, double);
    
    /* clean memory reserved for valist */
    va_end(valist);
    
    spline(num, ss, 0, dudt, 0.0, ts, uslog, b, c, d);
    uout = seval(num, t, ts, uslog, b, c, d);
    
    return(exp(uout));
}

/**
 * derivation of a spline function
 *
 * @param t point at which the spline should be evaluated
 * @param ti location of node i
 * @param pi spline value at node i
 * @param ss flag indicating whether slope at first node should be user defined
 * @param dudt user defined slope at first node
 *
 * @return dsplinedp(t)
 *
 */
double am_Dspline(int id, double t, int num, ...) {
    
    va_list valist;
    
    double uout;
    double ss;
    double dudt;
    
    double *ts = (double*) alloca(num*sizeof(double));
    double *us = (double*) alloca(num*sizeof(double));
    double *ps = (double*) alloca(num*sizeof(double));
    
    double *b = (double*) alloca(num*sizeof(double));
    double *c = (double*) alloca(num*sizeof(double));
    double *d = (double*) alloca(num*sizeof(double));
    
    int i;
    int j;
    int did;
    
    did = id/2 - 2;
    
    /* initialize valist for num number of arguments */
    va_start(valist, num);
    
    for (i=0; i<2*num; i+=2) {
        j = i/2;
        ts[j] = va_arg(valist, double);
        ps[j] = va_arg(valist, double);
        us[j] = 0.0;
    }
    us[did] = 1.0;
    ss = va_arg(valist, double);
    dudt = va_arg(valist, double);
    
    /* clean memory reserved for valist */
    va_end(valist);
    
    spline(num, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(num, t, ts, us, b, c, d);
    
    return(uout);
}

/**
 * derivation of an exponentiated spline function
 *
 * @param t point at which the spline should be evaluated
 * @param ti location of node i
 * @param pi spline value at node i
 * @param ss flag indicating whether slope at first node should be user defined
 * @param dudt user defined slope at first node
 *
 * @return dsplinedp(t)
 *
 */
double am_Dspline_pos(int id, double t, int num, ...) {
    
    va_list valist;
    
    double *ts = (double*) alloca(num*sizeof(double));
    double *us = (double*) alloca(num*sizeof(double));
    double *sus = (double*) alloca(num*sizeof(double));
    double *uslog = (double*) alloca(num*sizeof(double));

    double *b = (double*) alloca(num*sizeof(double));
    double *c = (double*) alloca(num*sizeof(double));
    double *d = (double*) alloca(num*sizeof(double));
    
    double uout;
    double ss;
    double dudt;
    double uspline_pos;
    double suspline;
    
    int i;
    int j;
    int did;
    
    did = id/2 - 2;
    
    /* initialize valist for num number of arguments */
    va_start(valist, num);
    
    for (i=0; i<2*num; i+=2) {
        j = i/2;
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
    
    spline(num, ss, 0, dudt, 0.0, ts, uslog, b, c, d);
    uspline_pos = exp(seval(num, t, ts, uslog, b, c, d));
    
    spline(num, ss, 0, dudt, 0.0, ts, sus, b, c, d);
    suspline = seval(num, t, ts, sus, b, c, d);
    uout = suspline * uspline_pos / us[did];
    
    return(uout);
}

/**
 * second derivation of a spline function
 *
 * @param t point at which the spline should be evaluated
 * @param ti location of node i
 * @param pi spline value at node i
 * @param ss flag indicating whether slope at first node should be user defined
 * @param dudt user defined slope at first node
 *
 * @return spline(t)
 *
 */
double am_DDspline(int id1, int id2, double t, int num, ...) {
    return(0.0);
}

/**
 * derivation of an exponentiated spline function
 *
 * @param t point at which the spline should be evaluated
 * @param ti location of node i
 * @param pi spline value at node i
 * @param ss flag indicating whether slope at first node should be user defined
 * @param dudt user defined slope at first node
 *
 * @return spline(t)
 *
 */
double am_DDspline_pos(int id1, int id2, double t, int num, ...) {
    
    va_list valist;
    
    double *ts = (double*) alloca(num*sizeof(double));
    double *us = (double*) alloca(num*sizeof(double));
    double *sus1 = (double*) alloca(num*sizeof(double));
    double *sus2 = (double*) alloca(num*sizeof(double));
    double *uslog = (double*) alloca(num*sizeof(double));
    
    double *b = (double*) alloca(num*sizeof(double));
    double *c = (double*) alloca(num*sizeof(double));
    double *d = (double*) alloca(num*sizeof(double));
    
    double uout;
    double ss;
    double dudt;
    double uspline_pos;
    double su1spline;
    double su2spline;
    
    int i;
    int j;
    int did1;
    int did2;
    
    did1 = id1/2 - 2;
    did2 = id2/2 - 2;
    
    /* initialize valist for num number of arguments */
    va_start(valist, num);
    
    for (i=0; i<2*num; i+=2) {
        j = i/2;
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
    
    spline(num, ss, 0, dudt, 0.0, ts, uslog, b, c, d);
    uspline_pos = exp(seval(num, t, ts, uslog, b, c, d));
    
    spline(num, ss, 0, dudt, 0.0, ts, sus1, b, c, d);
    su1spline = seval(num, t, ts, sus1, b, c, d);
    
    spline(num, ss, 0, dudt, 0.0, ts, sus2, b, c, d);
    su2spline = seval(num, t, ts, sus2, b, c, d);

    if (id1 == id2) {
        uout = (su1spline * su2spline - su1spline) * uspline_pos;
    }
    else {
        uout = su1spline * su2spline * uspline_pos;
    }
    uout = uout / us[did1] / us[did2];
    
    return(uout);
}
