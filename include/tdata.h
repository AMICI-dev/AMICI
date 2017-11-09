#ifndef AMICI_TDATA_H
#define AMICI_TDATA_H

#include <nvector/nvector_serial.h> /* defs. of serial NVECTOR fcts. and macros  */
#include <sundials/sundials_config.h>
#include <sundials/sundials_klu_impl.h> /* def. of type klu solver */
#include <sundials/sundials_math.h>     /* definition of ABS */
#include <sundials/sundials_sparse.h>   /* def. of type sparse stuff */
#include <sundials/sundials_types.h>    /* def. of type realtype */

namespace amici {

class UserData;
class ReturnData;
class Model;
class Solver;

/** @brief struct that provides temporary storage for different variables */
class TempData {

  public:
    TempData(const UserData *udata, Model *model, ReturnData *rdata);
    ~TempData();

    /**
     * @brief TempData is currently not copyable
     * @param other object to copy from
     */
    TempData (const TempData &other) = delete;

    /** parameter array, unscaled */
    realtype *p = nullptr;

    /** current time */
    realtype t;

    /** state vector */
    N_Vector x = nullptr;
    /** old state vector */
    N_Vector x_old = nullptr;
    /** array of state vectors at discontinuities*/
    N_Vector *x_disc = nullptr;
    /** array of differential state vectors at discontinuities*/
    N_Vector *xdot_disc = nullptr;
    /** array of old differential state vectors at discontinuities*/
    N_Vector *xdot_old_disc = nullptr;
    /** differential state vector */
    N_Vector dx = nullptr;
    /** old differential state vector */
    N_Vector dx_old = nullptr;
    /** time derivative state vector */
    N_Vector xdot = nullptr;
    /** old time derivative state vector */
    N_Vector xdot_old = nullptr;
    /** adjoint state vector */
    N_Vector xB = nullptr;
    /** old adjoint state vector */
    N_Vector xB_old = nullptr;
    /** differential adjoint state vector */
    N_Vector dxB = nullptr;
    /** quadrature state vector */
    N_Vector xQB = nullptr;
    /** old quadrature state vector */
    N_Vector xQB_old = nullptr;
    /** sensitivity state vector array */
    N_Vector *sx = nullptr;
    /** differential sensitivity state vector array */
    N_Vector *sdx = nullptr;
    /** Jacobian */
    DlsMat Jtmp = nullptr;

    /** parameter derivative of likelihood array */
    realtype *llhS0 = nullptr;
    /** data likelihood */
    realtype *Jy = nullptr;
    /** parameter derivative of data likelihood */
    realtype *dJydp = nullptr;
    /** observable derivative of data likelihood */
    realtype *dJydy = nullptr;
    /** observable sigma derivative of data likelihood */
    realtype *dJydsigma = nullptr;
    /** state derivative of data likelihood */
    realtype *dJydx = nullptr;
    /** event likelihood */
    realtype *Jz = nullptr;
    /** parameter derivative of event likelihood */
    realtype *dJzdp = nullptr;
    /** state derivative of event likelihood */
    realtype *dJzdx = nullptr;
    /** event ouput derivative of event likelihood */
    realtype *dJzdz = nullptr;
    /** event sigma derivative of event likelihood */
    realtype *dJzdsigma = nullptr;
    /** event ouput derivative of event likelihood at final timepoint */
    realtype *dJrzdz = nullptr;
    /** event sigma derivative of event likelihood at final timepoint */
    realtype *dJrzdsigma = nullptr;
    /** state derivative of event output */
    realtype *dzdx = nullptr;
    /** parameter derivative of event output */
    realtype *dzdp = nullptr;
    /** state derivative of event timepoint */
    realtype *drzdx = nullptr;
    /** parameter derivative of event timepoint */
    realtype *drzdp = nullptr;
    /** parameter derivative of observable */
    realtype *dydp = nullptr;
    /** state derivative of observable */
    realtype *dydx = nullptr;
    /** initial sensitivity of observable */
    realtype *yS0 = nullptr;
    /** data standard deviation */
    realtype *sigmay = nullptr;
    /** parameter derivative of data standard deviation */
    realtype *dsigmaydp = nullptr;
    /** event standard deviation */
    realtype *sigmaz = nullptr;
    /** parameter derivative of event standard deviation */
    realtype *dsigmazdp = nullptr;

    /** array of flags indicating which root has beend found.
     *  array of length nr with the indices of the user functions gi found to
     * have a
     *  root. For i = 0, . . . ,nr 1 if gi has a root, and = 0 if not.
     */
    int *rootsfound = nullptr;
    /** array of index which root has been found */
    int *rootidx = nullptr;
    /** array of number of found roots for a certain event type */
    int *nroots = nullptr;
    /** array of values of the root function */
    realtype *rootvals = nullptr;
    /** temporary rootval storage to check crossing in secondary event */
    realtype *rvaltmp = nullptr;

    /** flag indicating whether a certain heaviside function should be active or
       not */
    realtype *h = nullptr;

    /** change in x */
    realtype *deltax = nullptr;
    /** change in sx */
    realtype *deltasx = nullptr;
    /** change in xB */
    realtype *deltaxB = nullptr;
    /** change in qB */
    realtype *deltaqB = nullptr;

    /** integer for indexing of backwards problems */
    int which = 0;

    /** array containing the time-points of discontinuities*/
    realtype *discs = nullptr;
    /** array containing the index of discontinuities */
    realtype *irdiscs = nullptr;

    /** tempory storage of Jacobian data across functions */
    SlsMat J = nullptr;
    /** tempory storage of dxdotdp data across functions */
    realtype *dxdotdp = nullptr;
    /** tempory storage of w data across functions */
    realtype *w = nullptr;
    /** tempory storage of dwdx data across functions */
    realtype *dwdx = nullptr;
    /** tempory storage of dwdp data across functions */
    realtype *dwdp = nullptr;
    /** tempory storage of M data across functions */
    realtype *M = nullptr;
    /** tempory storage of dfdx data across functions */
    realtype *dfdx = nullptr;
    /** tempory storage of stau data across functions */
    realtype *stau = nullptr;

    /** number of parameters, copied from udata, necessary for deallocation */
    int nplist;

    /** current root index, will be increased during the forward solve and
     * decreased during backward solve */
    int iroot = 0;

    /** flag indicating whether a NaN in dxdotdp has been reported */
    booleantype nan_dxdotdp = false;
    /** flag indicating whether a NaN in J has been reported */
    booleantype nan_J = false;
    /** flag indicating whether a NaN in JDiag has been reported */
    booleantype nan_JDiag = false;
    /** flag indicating whether a NaN in JSparse has been reported */
    booleantype nan_JSparse = false;
    /** flag indicating whether a NaN in xdot has been reported */
    booleantype nan_xdot = false;
    /** flag indicating whether a NaN in xBdot has been reported */
    booleantype nan_xBdot = false;
    /** flag indicating whether a NaN in qBdot has been reported */
    booleantype nan_qBdot = false;

    /** attached UserData object */
    const UserData *udata = nullptr;
    /** attached Model object */
    Model *model = nullptr;
    /** attached ReturnData object */
    ReturnData *rdata = nullptr;
    /** attached Solver object */
    Solver *solver = nullptr;
};

} // namespace amici

#endif /* _MY_TDATA */
