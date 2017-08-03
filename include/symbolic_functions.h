#ifndef amici_symbolic_functions_h
#define amici_symbolic_functions_h

double amilog(double x);
double dirac(double x);
double heaviside(double x);

double am_min(double a, double b, double c);
double Dam_min(int id,double a, double b, double c);
double am_max(double a, double b, double c);
double Dam_max(int id,double a, double b, double c);

int amiIsNaN(double what);
int amiIsInf(double what);
double amiGetNaN();

/* sign */
double sign(double x);

/* splines */

double am_spline(double t, int num, ...);
double am_spline_pos(double t, int num, ...);
double am_Dspline(int id, double t, int num, ...);
double am_Dspline_pos(int id, double t, int num, ...);
double am_DDspline(int id1, int id2, double t, int num, ...);
double am_DDspline_pos(int id1, int id2, double t, int num, ...);

double spline3(double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt);
double spline_pos3(double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt);

double spline4(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt);
double spline_pos4(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt);

double spline5(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt);
double spline_pos5(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt);

double spline10(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt);
double spline_pos10(double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt);

double Dspline3(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt);
double Dspline_pos3(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt);

double Dspline4(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt);
double Dspline_pos4(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt);

double Dspline5(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt);
double Dspline_pos5(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt);

double Dspline10(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt);
double Dspline_pos10(int id, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt);

double DDspline3(int id1, int id2, double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt);
double DDspline_pos3(int id1, int id2, double t, double t1, double p1, double t2, double p2, double t3, double p3, int ss, double dudt);

double DDspline4(int id1, int id2, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt);
double DDspline_pos4(int id1, int id2, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, int ss, double dudt);

double DDspline5(int id1, int id2, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt);
double DDspline_pos5(int id1, int id2, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, int ss, double dudt);

double DDspline10(int id1, int id2, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt);
double DDspline_pos10(int id1, int id2, double t, double t1, double p1, double t2, double p2, double t3, double p3, double t4, double p4, double t5, double p5, double t6, double p6, double t7, double p7, double t8, double p8, double t9, double p9, double t10, double p10, int ss, double dudt);
#endif /* amici_symbolic_functions_h */
