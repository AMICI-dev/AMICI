#ifndef _MY_RDATA
#define _MY_RDATA
#include "udata.h"

/** @brief struct that stores all data which is later returned by the mex function
 *
 * NOTE: MATLAB stores multidimensional arrays in column-major order (FORTRAN-style)
 */
class ReturnData {
public:
    ReturnData();

    ReturnData(const UserData *udata);

    virtual ~ReturnData();

    /** timepoints */
    double *ts;
    /** time derivative */
    double *xdot;
    /** parameter derivative of time derivative */
    double *dxdotdp;
    /** state derivative of observables */
    double *dydx;
    /** parameter derivative of observables */
    double *dydp;
    /** Jacobian of differential equation right hand side */
    double *J;
    /** event output */
    double *z;
    /** event output sigma standard deviation */
    double *sigmaz;
    /** parameter derivative of event output */
    double *sz;
    /** parameter derivative of event output standard deviation */
    double *ssigmaz;
    /** event trigger output */
    double *rz;
    /** parameter derivative of event trigger output */
    double *srz;
    /** second order parameter derivative of event trigger output */
    double *s2rz;
    /** state */
    double *x;
    /** parameter derivative of state */
    double *sx;
    /** observable */
    double *y;
    /** observable standard deviation */
    double *sigmay;
    /** parameter derivative of observable */
    double *sy;
    /** parameter derivative of observable standard deviation */
    double *ssigmay;
    
    /** number of integration steps forward problem */
    double *numsteps;
    /** number of integration steps backward problem */
    double *numstepsS;
    /** number of right hand side evaluations forward problem */
    double *numrhsevals;
    /** number of right hand side evaluations backwad problem */
    double *numrhsevalsS;
    /** employed order forward problem */
    double *order;
    
    /** likelihood value */
    double *llh;
    /** chi2 value */
    double *chi2;
    /** parameter derivative of likelihood */
    double *sllh;
    /** second order parameter derivative of likelihood */
    double *s2llh;

    double *status;

protected:
    virtual void initFields(const UserData *udata);

    virtual void initField1(double **fieldPointer, const char *fieldName, int dim);

    /**
     * @ brief initialise matrix and attach to the field
     * @ param FIELD name of the field to which the matrix will be attached
     * @ param D1 number of rows in the matrix
     * @ param D2 number of columns in the matrix
     */

    virtual void initField2(double **fieldPointer, const char *fieldName, int dim1, int dim2);

    /**
     * @ brief initialise 3D tensor and attach to the field
     * @ param FIELD name of the field to which the tensor will be attached
     * @ param D1 number of rows in the tensor
     * @ param D2 number of columns in the tensor
     * @ param D3 number of elements in the third dimension of the tensor
     */

    virtual void initField3(double **fieldPointer, const char *fieldName, int dim1, int dim2, int dim3);

    /**
     * @ brief initialise 4D tensor and attach to the field
     * @ param FIELD name of the field to which the tensor will be attached
     * @ param D1 number of rows in the tensor
     * @ param D2 number of columns in the tensor
     * @ param D3 number of elements in the third dimension of the tensor
     * @ param D4 number of elements in the fourth dimension of the tensor
     */

    virtual void initField4(double **fieldPointer, const char *fieldName, int dim1, int dim2, int dim3, int dim4);

    
};

#endif /* _MY_RDATA */
