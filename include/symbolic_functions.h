#ifndef amici_symbolic_functions_h
#define amici_symbolic_functions_h
#include <math.h>
#include <include/spline.h>
#include <src/symbolic_functions.c>

/* ge */

static int cw_ge(double a, double b);
static int Dcw_ge(int id, double a, double b);

/* gt */

static int cw_gt(double a, double b);

static int Dcw_gt(int id, double a, double b);

/* le */

static int cw_le(double a, double b);
static int Dcw_le(int id, double a, double b);

/* lt */

static int cw_lt(double a, double b);

static int Dcw_lt(int id, double a, double b);

/* If */

static double cw_if(int condition, double truepart, double falsepart);
static double Dcw_if(int id, int condition, double truepart, double falsepart);

/* min */

static double cw_min(double a, double b);
static double Dcw_min(int id,double a, double b);

static double cw_max(double a, double b);
static double Dcw_max(int id,double a, double b);

/* and */

static int cw_and(int a, int b);
static double Dcw_and(int id,int a, int b);

/* or */

static int cw_or(int a, int b);
static double Dcw_or(int id,int a, int b);


/* step */

static double heaviside(double x);
static double dirac(double x);
static double Ddirac(int ip, double x);

/* sign */
static double step(double x);

/* splines */

static double spline3(double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt);
static double spline_pos3(double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt);

static double spline4(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt);
static double spline_pos4(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt);

static double spline5(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt);
static double spline_pos5(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt);

static double spline10(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt);
static double spline_pos10(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt);

static double Dspline3(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt);
static double Dspline_pos3(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt);

static double Dspline4(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt);
static double Dspline_pos4(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt);

static double Dspline5(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt);
static double Dspline_pos5(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt);

static double Dspline10(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt);
static double Dspline_pos10(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt);
#endif /* amici_symbolic_functions_h */