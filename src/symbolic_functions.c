/**
 * @file   symbolic_functions.c
 * @brief  definition of symbolic functions
 *
 * This file contains definitions of various symbolic functions which
 */


#if (_MSC_VER >= 1000)
#define fmax(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
#define fmin(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })
#endif

#include <math.h>
#include <mex.h>
#include <float.h>
#include <include/spline.h>

#undef ts

/*! bool return value true */
#define TRUE 1
/*! bool return value false */
#define FALSE 0


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
 * @return if(a < b) then a else b @type double
 *
 */
double am_min(double a, double b) {
    return(fmin(a,b));
}

/**
 * parameter derivative of c implementation of matlab function min
 *
 * @param id argument index for differentiation
 * @param a bool1 @type double
 * @param b bool2 @type double
 * @return id == 1:  if(a < b) then 1 else 0 @type double
 * @return id == 2:  if(a < b) then 0 else 1 @type double
 *
 */
double Dam_min(int id,double a, double b) {
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
 * @return if(a > b) then a else b @type double
 *
 */
double am_max(double a, double b) {
    return(fmax(a,b));
}

/**
 * parameter derivative of c implementation of matlab function max
 *
 * @param id argument index for differentiation
 * @param a bool1 @type double
 * @param b bool2 @type double
 * @return id == 1:  if(a > b) then 1 else 0 @type double
 * @return id == 2:  if(a > b) then 0 else 1 @type double
 *
 */
double Dam_max(int id,double a, double b) {
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
 * spline function with 3 nodes
 *
 * @param t point at which the spline should be evaluated
 * @param t1 location of node 1
 * @param p1 spline value at node 1
 * @param t2 location of node 2
 * @param p2 spline value at node 2
 * @param t3 location of node 3
 * @param p3 spline value at node 3
 * @param ss flag indicating whether slope at first node should be user defined
 * @param dudt user defined slope at first node
 *
 * @return spline(t)
 *
 */
double spline3(double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt) {
    double uout;
    
    double ts[3];
    double us[3];
    
    double b[3];
    double c[3];
    double d[3];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    
    spline(3, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(3, t, ts, us, b, c, d);
    
    return(uout);
}


/**
 * positive spline function with 3 nodes
 *
 * @param t point at which the spline should be evaluated
 * @param t1 location of node 1
 * @param p1 spline value at node 1
 * @param t2 location of node 2
 * @param p2 spline value at node 2
 * @param t3 location of node 3
 * @param p3 spline value at node 3
 * @param ss flag indicating whether slope at first node should be user defined
 * @param dudt user defined slope at first node
 *
 * @return spline(t)
 *
 */
double spline_pos3(double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt) {
    int is;
    double uout;
    
    double ts[3];
    double us[3];
    double uslog[3];
    
    double b[3];
    double c[3];
    double d[3];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    
    for (is = 0; is<3; is++){
        uslog[is] = log(us[is]);
    }
    
    spline(3, ss, 0, dudt, 0.0, ts, uslog, b, c, d);
    uout = seval(3, t, ts, uslog, b, c, d);
    
    return(exp(uout));
}

/**
 * spline function with 4 nodes
 *
 * @param t point at which the spline should be evaluated
 * @param t1 location of node 1
 * @param p1 spline value at node 1
 * @param t2 location of node 2
 * @param p2 spline value at node 2
 * @param t3 location of node 3
 * @param p3 spline value at node 3
 * @param t4 location of node 4
 * @param p4 spline value at node 4
 * @param ss flag indicating whether slope at first node should be user defined
 * @param dudt user defined slope at first node
 *
 * @return spline(t)
 *
 */
double spline4(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt) {
    double uout;
    
    double ts[4];
    double us[4];
    
    double b[4];
    double c[4];
    double d[4];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    
    spline(4, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(4, t, ts, us, b, c, d);
    
    return(uout);
}

/**
 * positive spline function with 4 nodes
 *
 * @param t point at which the spline should be evaluated
 * @param t1 location of node 1
 * @param p1 spline value at node 1
 * @param t2 location of node 2
 * @param p2 spline value at node 2
 * @param t3 location of node 3
 * @param p3 spline value at node 3
 * @param t4 location of node 4
 * @param p4 spline value at node 4
 * @param ss flag indicating whether slope at first node should be user defined
 * @param dudt user defined slope at first node
 *
 * @return spline(t)
 *
 */
double spline_pos4(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt) {
    int is;
    double uout;
    
    double ts[4];
    double us[4];
    double uslog[4];
    
    double b[4];
    double c[4];
    double d[4];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    
    for (is = 0; is<4; is++){
        uslog[is] = log(us[is]);
    }
    
    spline(4, ss, 0, dudt, 0.0, ts, uslog, b, c, d);
    uout = seval(4, t, ts, uslog, b, c, d);
    
    return(exp(uout));
}

/**
 * spline function with 5 nodes
 *
 * @param t point at which the spline should be evaluated
 * @param t1 location of node 1
 * @param p1 spline value at node 1
 * @param t2 location of node 2
 * @param p2 spline value at node 2
 * @param t3 location of node 3
 * @param p3 spline value at node 3
 * @param t4 location of node 4
 * @param p4 spline value at node 4
 * @param t5 location of node 5
 * @param p5 spline value at node 5
 * @param ss flag indicating whether slope at first node should be user defined
 * @param dudt user defined slope at first node
 *
 * @return spline(t)
 *
 */
double spline5(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt) {
    double uout;
    
    double ts[5];
    double us[5];
    
    double b[5];
    double c[5];
    double d[5];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    us[4] = p5;
    
    spline(5, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(5, t, ts, us, b, c, d);
    
    return(uout);
}

/**
 * positive spline function with 5 nodes
 *
 * @param t point at which the spline should be evaluated
 * @param t1 location of node 1
 * @param p1 spline value at node 1
 * @param t2 location of node 2
 * @param p2 spline value at node 2
 * @param t3 location of node 3
 * @param p3 spline value at node 3
 * @param t4 location of node 4
 * @param p4 spline value at node 4
 * @param t5 location of node 5
 * @param p5 spline value at node 5
 * @param ss flag indicating whether slope at first node should be user defined
 * @param dudt user defined slope at first node
 *
 * @return spline(t)
 *
 */
double spline_pos5(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt) {
    int is;
    double uout;
    
    double ts[5];
    double us[5];
    double uslog[5];
    
    double b[5];
    double c[5];
    double d[5];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    us[4] = p5;
    
    for (is = 0; is<5; is++){
        uslog[is] = log(us[is]);
    }
    
    spline(5, ss, 0, dudt, 0.0, ts, uslog, b, c, d);
    uout = seval(5, t, ts, uslog, b, c, d);
    
    return(exp(uout));
}

/**
 * spline function with 10 nodes
 *
 * @param t point at which the spline should be evaluated
 * @param t1 location of node 1
 * @param p1 spline value at node 1
 * @param t2 location of node 2
 * @param p2 spline value at node 2
 * @param t3 location of node 3
 * @param p3 spline value at node 3
 * @param t4 location of node 4
 * @param p4 spline value at node 4
 * @param t5 location of node 5
 * @param p5 spline value at node 5
 * @param t6 location of node 6
 * @param p6 spline value at node 6
 * @param t7 location of node 7
 * @param p7 spline value at node 7
 * @param t8 location of node 8
 * @param p8 spline value at node 8
 * @param t9 location of node 9
 * @param p9 spline value at node 9
 * @param t10 location of node 10
 * @param p10 spline value at node 10
 * @param ss flag indicating whether slope at first node should be user defined
 * @param dudt user defined slope at first node
 *
 * @return spline(t)
 *
 */
double spline10(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt) {
    double uout;
    
    double ts[10];
    double us[10];
    
    double b[10];
    double c[10];
    double d[10];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    ts[5] = t6;
    ts[6] = t7;
    ts[7] = t8;
    ts[8] = t9;
    ts[9] = t10;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    us[4] = p5;
    us[5] = p6;
    us[6] = p7;
    us[7] = p8;
    us[8] = p9;
    us[9] = p10;
    
    spline(10, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(10, t, ts, us, b, c, d);
    
    return(uout);
}

/**
 * positive spline function with 10 nodes
 *
 * @param t point at which the spline should be evaluated
 * @param t1 location of node 1
 * @param p1 spline value at node 1
 * @param t2 location of node 2
 * @param p2 spline value at node 2
 * @param t3 location of node 3
 * @param p3 spline value at node 3
 * @param t4 location of node 4
 * @param p4 spline value at node 4
 * @param t5 location of node 5
 * @param p5 spline value at node 5
 * @param t6 location of node 6
 * @param p6 spline value at node 6
 * @param t7 location of node 7
 * @param p7 spline value at node 7
 * @param t8 location of node 8
 * @param p8 spline value at node 8
 * @param t9 location of node 9
 * @param p9 spline value at node 9
 * @param t10 location of node 10
 * @param p10 spline value at node 10
 * @param ss flag indicating whether slope at first node should be user defined
 * @param dudt user defined slope at first node
 *
 * @return spline(t)
 *
 */
double spline_pos10(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt) {
    int is;
    double uout;
    
    double ts[10];
    double us[10];
    double uslog[10];
    
    double b[10];
    double c[10];
    double d[10];
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    ts[5] = t6;
    ts[6] = t7;
    ts[7] = t8;
    ts[8] = t9;
    ts[9] = t10;
    
    us[0] = p1;
    us[1] = p2;
    us[2] = p3;
    us[3] = p4;
    us[4] = p5;
    us[5] = p6;
    us[6] = p7;
    us[7] = p8;
    us[8] = p9;
    us[9] = p10;
    
    for (is = 0; is<10; is++){
        uslog[is] = log(us[is]);
    }
    
    spline(10, ss, 0, dudt, 0.0, ts, uslog, b, c, d);
    uout = seval(10, t, ts, uslog, b, c, d);
    
    return(exp(uout));
}

/**
 * parameter derivative of spline function with 3 nodes
 *
 * @param id argument index for differentiation
 * @param t point at which the spline should be evaluated
 * @param t1 location of node 1
 * @param p1 spline value at node 1
 * @param t2 location of node 2
 * @param p2 spline value at node 2
 * @param t3 location of node 3
 * @param p3 spline value at node 3
 * @param ss flag indicating whether slope at first node should be user defined
 * @param dudt user defined slope at first node
 *
 * @return dspline(t)dp(id)
 *
 */
double Dspline3(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt) {
    double uout;
    
    double ts[3];
    double us[3];
    
    double b[3];
    double c[3];
    double d[3];
    int did;
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    
    us[0] = 0.0;
    us[1] = 0.0;
    us[2] = 0.0;
    
    did = floor(id/2-1);
    us[did] = 1.0;
    
    spline(3, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(3, t, ts, us, b, c, d);
    
    return(uout);
}

/**
 * parameter derivative of positive spline function with 3 nodes
 *
 * @param id argument index for differentiation
 * @param t point at which the spline should be evaluated
 * @param t1 location of node 1
 * @param p1 spline value at node 1
 * @param t2 location of node 2
 * @param p2 spline value at node 2
 * @param t3 location of node 3
 * @param p3 spline value at node 3
 * @param ss flag indicating whether slope at first node should be user defined
 * @param dudt user defined slope at first node
 *
 * @return dspline(t)dp(id)
 *
 */
double Dspline_pos3(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt) {
    
    double uout;
    double uspline_pos3;
    double suspline3;
    
    double ps[3];
    int did;
    
    
    ps[0] = p1;
    ps[1] = p2;
    ps[2] = p3;
    
    did = floor(id/2-1);
    
    uspline_pos3 = spline_pos3(t, t1, p1, t2, p2, t3, p3, ss, dudt);
    suspline3 = Dspline3(id, t, t1, p1, t2, p2, t3, p3, ss, dudt);
    uout = suspline3 * uspline_pos3 * log(10.0);
    uout = uout / ps[did] / log(10.0);
    
    return(uout);
}

/**
 * parameter derivative of spline function with 4 nodes
 *
 * @param id argument index for differentiation
 * @param t point at which the spline should be evaluated
 * @param t1 location of node 1
 * @param p1 spline value at node 1
 * @param t2 location of node 2
 * @param p2 spline value at node 2
 * @param t3 location of node 3
 * @param p3 spline value at node 3
 * @param t4 location of node 4
 * @param p4 spline value at node 4
 * @param ss flag indicating whether slope at first node should be user defined
 * @param dudt user defined slope at first node
 *
 * @return dspline(t)dp(id)
 *
 */
double Dspline4(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt) {
    double uout;
    
    double ts[4];
    double us[4];
    
    double b[4];
    double c[4];
    double d[4];
    int did;
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    
    us[0] = 0.0;
    us[1] = 0.0;
    us[2] = 0.0;
    us[3] = 0.0;
    
    did = floor(id/2-1);
    us[did] = 1.0;
    
    spline(4, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(4, t, ts, us, b, c, d);
    
    return(uout);
}

/**
 * parameter derivative of positive spline function with 4 nodes
 *
 * @param id argument index for differentiation
 * @param t point at which the spline should be evaluated
 * @param t1 location of node 1
 * @param p1 spline value at node 1
 * @param t2 location of node 2
 * @param p2 spline value at node 2
 * @param t3 location of node 3
 * @param p3 spline value at node 3
 * @param t4 location of node 4
 * @param p4 spline value at node 4
 * @param ss flag indicating whether slope at first node should be user defined
 * @param dudt user defined slope at first node
 *
 * @return dspline(t)dp(id)
 *
 */
double Dspline_pos4(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt) {
    
    double uout;
    double uspline_pos4;
    double suspline4;
    
    double ps[4];
    int did;
    
    ps[0] = p1;
    ps[1] = p2;
    ps[2] = p3;
    ps[3] = p4;
    
    did = floor(id/2-1);
    
    uspline_pos4 = spline_pos4(t, t1, p1, t2, p2, t3, p3, t4, p4, ss, dudt);
    suspline4 = Dspline4(id, t, t1, p1, t2, p2, t3, p3, t4, p4, ss, dudt);
    uout = suspline4 * uspline_pos4 * log(10.0);
    uout = uout / ps[did] / log(10.0);
    
    return(uout);
}

/**
 * parameter derivative of spline function with 5 nodes
 *
 * @param id argument index for differentiation
 * @param t point at which the spline should be evaluated
 * @param t1 location of node 1
 * @param p1 spline value at node 1
 * @param t2 location of node 2
 * @param p2 spline value at node 2
 * @param t3 location of node 3
 * @param p3 spline value at node 3
 * @param t4 location of node 4
 * @param p4 spline value at node 4
 * @param t5 location of node 5
 * @param p5 spline value at node 5
 * @param ss flag indicating whether slope at first node should be user defined
 * @param dudt user defined slope at first node
 *
 * @return dspline(t)dp(id)
 *
 */
double Dspline5(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt) {
    double uout;
    
    double ts[5];
    double us[5];
    
    double b[5];
    double c[5];
    double d[5];
    int did;
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    
    us[0] = 0.0;
    us[1] = 0.0;
    us[2] = 0.0;
    us[3] = 0.0;
    us[4] = 0.0;
    
    did = floor(id/2-1);
    us[did] = 1.0;
    
    spline(5, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(5, t, ts, us, b, c, d);
    
    return(uout);
}

/**
 * parameter derivative of positive spline function with 5 nodes
 *
 * @param id argument index for differentiation
 * @param t point at which the spline should be evaluated
 * @param t1 location of node 1
 * @param p1 spline value at node 1
 * @param t2 location of node 2
 * @param p2 spline value at node 2
 * @param t3 location of node 3
 * @param p3 spline value at node 3
 * @param t4 location of node 4
 * @param p4 spline value at node 4
 * @param t5 location of node 5
 * @param p5 spline value at node 5
 * @param ss flag indicating whether slope at first node should be user defined
 * @param dudt user defined slope at first node
 *
 * @return dspline(t)dp(id)
 *
 */
double Dspline_pos5(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt) {
    
    double uout;
    double uspline_pos5;
    double suspline5;
    
    double ps[5];
    int did;
    
    ps[0] = p1;
    ps[1] = p2;
    ps[2] = p3;
    ps[3] = p4;
    ps[4] = p5;
    
    did = floor(id/2-1);
    
    uspline_pos5 = spline_pos5(t, t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, ss, dudt);
    suspline5 = Dspline5(id, t, t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, ss, dudt);
    uout = suspline5 * uspline_pos5 * log(10.0);
    uout = uout / ps[did] / log(10.0);
    
    return(uout);
}

/**
 * parameter derivative of spline function with 10 nodes
 *
 * @param id argument index for differentiation
 * @param t point at which the spline should be evaluated
 * @param t1 location of node 1
 * @param p1 spline value at node 1
 * @param t2 location of node 2
 * @param p2 spline value at node 2
 * @param t3 location of node 3
 * @param p3 spline value at node 3
 * @param t4 location of node 4
 * @param p4 spline value at node 4
 * @param t5 location of node 5
 * @param p5 spline value at node 5
 * @param t6 location of node 6
 * @param p6 spline value at node 6
 * @param t7 location of node 7
 * @param p7 spline value at node 7
 * @param t8 location of node 8
 * @param p8 spline value at node 8
 * @param t9 location of node 9
 * @param p9 spline value at node 9
 * @param t10 location of node 10
 * @param p10 spline value at node 10
 * @param ss flag indicating whether slope at first node should be user defined
 * @param dudt user defined slope at first node
 *
 * @return dspline(t)dp(id)
 *
 */
double Dspline10(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt) {
    double uout;
    
    double ts[10];
    double us[10];
    
    double b[10];
    double c[10];
    double d[10];
    int did;
    
    ts[0] = t1;
    ts[1] = t2;
    ts[2] = t3;
    ts[3] = t4;
    ts[4] = t5;
    ts[5] = t6;
    ts[6] = t7;
    ts[7] = t8;
    ts[8] = t9;
    ts[9] = t10;
    
    us[0] = 0.0;
    us[1] = 0.0;
    us[2] = 0.0;
    us[3] = 0.0;
    us[4] = 0.0;
    us[5] = 0.0;
    us[6] = 0.0;
    us[7] = 0.0;
    us[8] = 0.0;
    us[9] = 0.0;
    
    did = floor(id/2-1);
    us[did] = 1.0;
    
    spline(10, ss, 0, dudt, 0.0, ts, us, b, c, d);
    uout = seval(10, t, ts, us, b, c, d);
    
    return(uout);
}

/**
 * parameter derivative of positive spline function with 10 nodes
 *
 * @param id argument index for differentiation
 * @param t point at which the spline should be evaluated
 * @param t1 location of node 1
 * @param p1 spline value at node 1
 * @param t2 location of node 2
 * @param p2 spline value at node 2
 * @param t3 location of node 3
 * @param p3 spline value at node 3
 * @param t4 location of node 4
 * @param p4 spline value at node 4
 * @param t5 location of node 5
 * @param p5 spline value at node 5
 * @param t6 location of node 6
 * @param p6 spline value at node 6
 * @param t7 location of node 7
 * @param p7 spline value at node 7
 * @param t8 location of node 8
 * @param p8 spline value at node 8
 * @param t9 location of node 9
 * @param p9 spline value at node 9
 * @param t10 location of node 10
 * @param p10 spline value at node 10
 * @param ss flag indicating whether slope at first node should be user defined
 * @param dudt user defined slope at first node
 *
 * @return dspline(t)dp(id)
 *
 */
double Dspline_pos10(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt) {
    
    double uout;
    double uspline_pos10;
    double suspline10;
    
    double ps[10];
    int did;
    
    ps[0] = p1;
    ps[1] = p2;
    ps[2] = p3;
    ps[3] = p4;
    ps[4] = p5;
    ps[5] = p6;
    ps[6] = p7;
    ps[7] = p8;
    ps[8] = p9;
    ps[9] = p10;
    
    did = floor(id/2-1);
    
    uspline_pos10 = spline_pos10(t, t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6, t7, p7, t8, p8, t9, p9, t10, p10, ss, dudt);
    suspline10 = Dspline10(id, t, t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6, t7, p7, t8, p8, t9, p9, t10, p10, ss, dudt);
    uout = suspline10 * uspline_pos10 * log(10.0);
    uout = uout / ps[did] / log(10.0);
    
    return(uout);
}