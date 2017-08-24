#ifndef MODEL_H
#define MODEL_H

#include <include/amici.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_direct.h>
#include <sundials/sundials_sparse.h>

class UserData;
class ExpData;

/**
 * @brief The Model class represents an AMICI ODE model.
 * The model does not contain any data, its state should never change.
 */
class Model {
  public:
    Model()
        : np(0), nk(0), nx(0), nxtrue(0), ny(0), nytrue(0), nz(0), nztrue(0),
          ne(0), nw(0), ndwdx(0), ndwdp(0), nnz(0), nJ(0), ubw(0), lbw(0),
          o2mode(AMICI_O2MODE_NONE) {}

    Model(int np, int nx, int nxtrue, int nk, int ny, int nytrue, int nz,
          int nztrue, int ne, int nJ, int nw, int ndwdx, int ndwdp, int nnz,
          int ubw, int lbw, AMICI_o2mode o2mode)
        : np(np), nk(nk), nx(nx), nxtrue(nxtrue), ny(ny), nytrue(nytrue),
          nz(nz), nztrue(nztrue), ne(ne), nw(nw), ndwdx(ndwdx), ndwdp(ndwdp),
          nnz(nnz), nJ(nJ), ubw(ubw), lbw(lbw), o2mode(o2mode) {}

    virtual Solver *getSolver() = 0;

    /** Initial states **/
    virtual int fx0(N_Vector x0, void *user_data) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Initial value for time derivative of states **/
    virtual int fdx0(N_Vector x0, N_Vector dx0, void *user_data) { return AMICI_SUCCESS; }

//    virtual int fdx0(const realtype *k, realtype *x) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Initial value for time derivative of states **/
    virtual int fsx0(N_Vector *sx0, N_Vector x, N_Vector dx, void *user_data) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Sensitivity of initial states x0 */
    virtual int fsdx0(N_Vector *sdx0, N_Vector x, N_Vector dx, void *user_data) { return AMICI_SUCCESS; }

    /** Jacobian of xdot with respect to states x */
    virtual int fJ(long int N, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Jacobian of xBdot with respect to adjoint state xB */
    virtual int fJB(long int NeqBdot, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** diagonalized Jacobian (for preconditioning) */
    virtual int fJDiag(realtype t, N_Vector JDiag, N_Vector x, void *user_data)  { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Matrix vector product of J with a vector v (for iterative solvers) */
    virtual int fJv(N_Vector v, N_Vector Jv, realtype t, N_Vector x, N_Vector xdot, void *user_data, N_Vector tmp)  { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Event trigger function for events */
    virtual int froot(realtype t, N_Vector x, N_Vector dx, realtype *root, void *user_data) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Event root function of events (equal to froot but does not include non-output events) */
    virtual int frz(realtype t, int ie, N_Vector x, TempData *tdata, ReturnData *rdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Sensitivity of rz, total derivative */
    virtual int fsrz(realtype t, int ie, N_Vector x, N_Vector *sx, TempData *tdata, ReturnData *rdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Sensitivity of event timepoint, total derivative */
    virtual int fstau(realtype t, int ie, N_Vector x, N_Vector *sx, TempData *tdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Observables / measurements */
    virtual int fy(realtype t, int it, N_Vector x, void *user_data, ReturnData *rdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Sensitivity of observables y w.r.t. model parameters p */
    virtual int fdydp(realtype t, int it, N_Vector x, TempData *tdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Sensitivity of observables y w.r.t. state variables x */
    virtual int fdydx(realtype t, int it, N_Vector x, TempData *tdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Event-resolved measurements */
    virtual int fz(realtype t, int ie, N_Vector x, TempData *tdata, ReturnData *rdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Sensitivity of z, total derivative */
    virtual int fsz(realtype t, int ie, N_Vector x, N_Vector *sx, TempData *tdata, ReturnData *rdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Sensitivity of event-resolved measurements z w.r.t. to model parameters p */
    virtual int fdzdp(realtype t, int ie, N_Vector x, TempData *tdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Sensitivity of event-resolved measurements z w.r.t. to model states x */
    virtual int fdzdx(realtype t, int ie, N_Vector x, TempData *tdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Sensitivity of event-resolved measurements rz w.r.t. to model parameters p */
    virtual int fdrzdp(realtype t, int ie, N_Vector x, TempData *tdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Sensitivity of event-resolved measurements rz w.r.t. to model states x */
    virtual int fdrzdx(realtype t, int ie, N_Vector x, TempData *tdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Right hand side of differential equation for states x */
    virtual int fxdot(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, void *user_data) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Right hand side of differential equation for adjoint state xB */
    virtual int fxBdot(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Right hand side of integral equation for quadrature states qB */
    virtual int fqBdot(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Sensitivity of dx/dt w.r.t. model parameters p */
    virtual int fdxdotdp(realtype t, N_Vector x, N_Vector dx, void *user_data) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** State update functions for events */
    virtual int fdeltax(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, TempData *tdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Sensitivity update functions for events, total derivative */
    virtual int fdeltasx(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, N_Vector *sx, TempData *tdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Adjoint state update functions for events */
    virtual int fdeltaxB(realtype t, int ie, N_Vector x, N_Vector xB, N_Vector xdot, N_Vector xdot_old, TempData *tdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Quadrature state update functions for events */
    virtual int fdeltaqB(realtype t, int ie, N_Vector x, N_Vector xB, N_Vector qBdot, N_Vector xdot, N_Vector xdot_old, TempData *tdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Standard deviation of measurements */
    virtual int fsigma_y(realtype t, TempData *tdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Sensitivity of standard deviation of measurements w.r.t. model parameters p */
    virtual int fdsigma_ydp(realtype t, TempData *tdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Standard deviation of events */
    virtual int fsigma_z(realtype t, int ie, TempData *tdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Sensitivity of standard deviation of events w.r.t. model parameters p */
    virtual int fdsigma_zdp(realtype t, int ie, TempData *tdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** negative log-likelihood of time-resolved measurements y */
    virtual int fJy(realtype t, int it, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** negative log-likelihood of event-resolved measurements z */
    virtual int fJz(realtype t, int ie, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** regularization of negative log-likelihood with roots of event-resolved measurements rz */
    virtual int fJrz(realtype t, int ie, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Sensitivity of time-resolved measurement negative log-likelihood Jy w.r.t. observables y */
    virtual int fdJydy(realtype t, int it, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Sensitivity of time-resolved measurement negative log-likelihood Jy w.r.t. standard deviation sigma */
    virtual int fdJydsigma(realtype t, int it, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Sensitivity of event-resolved measurement negative log-likelihood Jz w.r.t. event observables z */
    virtual int fdJzdz(realtype t, int ie, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Sensitivity of event-resolved measurement negative log-likelihood Jz w.r.t. standard deviation sigma */
    virtual int fdJzdsigma(realtype t, int ie, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Sensitivity of event-resolved measurement negative log-likelihood regularization Jrz w.r.t. event observables z */
    virtual int fdJrzdz(realtype t, int ie, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Sensitivity of event-resolved measurement negative log-likelihood regularization Jrz w.r.t. standard deviation sigma */
    virtual int fdJrzdsigma(realtype t, int ie, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Right hand side of differential equation for state sensitivities sx */
    virtual int fsxdot(int Ns, realtype t, N_Vector x, N_Vector xdot,int ip,  N_Vector sx, N_Vector sxdot, void *user_data, N_Vector tmp1, N_Vector tmp2) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** J in sparse form (for sparse solvers from the SuiteSparse Package) */
    virtual int fJSparse(realtype t, N_Vector x, N_Vector xdot, SlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** J in banded form (for banded solvers) */
    virtual int fJBand(long int N, long int mupper, long int mlower, realtype t, N_Vector x, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** JB in banded form (for banded solvers) */
    virtual int fJBandB(long int NeqBdot, long int mupper, long int mlower, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** Matrix vector product of JB with a vector v (for iterative solvers) */
    virtual int fJvB(N_Vector vB, N_Vector JvB, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data, N_Vector tmpB) { return AMICI_ERROR_NOT_IMPLEMENTED; }

    /** JB in sparse form (for sparse solvers from the SuiteSparse Package) */
    virtual int fJSparseB(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, SlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) { return AMICI_ERROR_NOT_IMPLEMENTED; }


    virtual ~Model();

    // Generic implementations
    /** Sensitivity of measurements y, total derivative */
    int fsy(int it, TempData *tdata, ReturnData *rdata);

    /** Sensitivity of z at final timepoint (ignores sensitivity of timepoint), total derivative *//
    int fsz_tf(int ie, TempData *tdata, ReturnData *rdata);

    /* Sensitivity of time-resolved measurement negative log-likelihood Jy, total derivative */
    int fsJy(int it, TempData *tdata, ReturnData *rdata);

    /* Sensitivity of time-resolved measurement negative log-likelihood Jy w.r.t. parameters  */
    int fdJydp(int it, TempData *tdata, const ExpData *edata, ReturnData *rdata);

    /* Sensitivity of time-resolved measurement negative log-likelihood Jy w.r.t. state variables  */
    int fdJydx(int it, TempData *tdata, const ExpData *edata);

    /* Sensitivity of event-resolved measurement negative log-likelihood Jz, total derivative */
    int fsJz(int ie, TempData *tdata, const ReturnData *rdata);

    /* Sensitivity of event-resolved measurement negative log-likelihood Jz w.r.t. parameters  */
    int fdJzdp(int ie, TempData *tdata, const ExpData *edata, ReturnData *rdata);

    /* Sensitivity of event-resolved measurement negative log-likelihood Jz w.r.t. state variables  */
    int fdJzdx(int ie, TempData *tdata, const ExpData *edata);

    /* initialization of model properties  */
    int initialize(UserData *udata, TempData *tdata);

    /* initialization of initial states  */
    int initializeStates(double *x0data, TempData *tdata);

    /* initialization of heaviside function variables  */
    int initHeaviside(TempData *tdata);

    /* Model dimensions */
    /** total number of model parameters */
    const int np;
    /** number of fixed parameters */
    const int nk;
    /** number of states */
    const int nx;
    /** number of states in the unaugmented system */
    const int nxtrue;
    /** number of observables */
    const int ny;
    /** number of observables in the unaugmented system */
    const int nytrue;
    /** number of event outputs */
    const int nz;
    /** number of event outputs in the unaugmented system */
    const int nztrue;
    /** number of events */
    const int ne;
    /** number of common expressions */
    const int nw;
    /** number of derivatives of common expressions wrt x */
    const int ndwdx;
    /** number of derivatives of common expressions wrt p */
    const int ndwdp;
    /** number of nonzero entries in jacobian */
    const int nnz;
    /** dimension of the augmented objective function for 2nd order ASA */
    const int nJ;
    /** upper bandwith of the jacobian */
    const int ubw;
    /** lower bandwith of the jacobian */
    const int lbw;
    /** flag indicating whether for sensi == AMICI_SENSI_ORDER_SECOND
     * directional or full second order derivative will be computed */
    const AMICI_o2mode o2mode;
    /** index indicating to which event an event output belongs */
    int *z2event = nullptr;
    /** flag array for DAE equations */
    realtype *idlist = nullptr;
};

#endif // MODEL_H
