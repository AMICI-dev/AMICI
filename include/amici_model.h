#ifndef AMICI_fH
#define AMICI_fH

#include <include/amici_exception.h>
#include <include/amici_defines.h>
#include <include/amici_vector.h>
#include <include/symbolic_functions.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_direct.h>
#include <sundials/sundials_sparse.h>

#include <vector>

namespace amici {
    
    class UserData;
    class ReturnData;
    class ExpData;
    class Solver;

    /**
     * @brief The Model class represents an AMICI ODE model.
     * The model can compute various model related quantities based
     * on symbolically generated code.
     */
    class Model {
    public:
        /** default constructor */
        Model()
        : nx(0), nxtrue(0), ny(0), nytrue(0), nz(0), nztrue(0),
        ne(0), nw(0), ndwdx(0), ndwdp(0), nnz(0), nJ(0), ubw(0), lbw(0),
        o2mode(AMICI_O2MODE_NONE) {}
        
        /** constructor with model dimensions
         * @param nx number of state variables
         * @param nxtrue number of state variables of the non-augmented model
         * @param ny number of observables
         * @param nytrue number of observables of the non-augmented model
         * @param nz number of event observables
         * @param nztrue number of event observables of the non-augmented model
         * @param ne number of events
         * @param nJ number of objective functions
         * @param nw number of repeating elements
         * @param ndwdx number of nonzero elements in the x derivative of the
         * repeating elements
         * @param ndwdp number of nonzero elements in the p derivative of the
         * repeating elements
         * @param nnz number of nonzero elements in Jacobian
         * @param ubw upper matrix bandwidth in the Jacobian
         * @param lbw lower matrix bandwidth in the Jacobian
         * @param o2mode second order sensitivity mode
         * @param p parameters
         * @param k constants
         * @param plist indexes wrt to which sensitivities are to be computed
         * @param idlist indexes indicating algebraic components (DAE only)
         * @param z2event mapping of event outputs to events
         */
        Model(const int nx, const int nxtrue,
              const int ny, const int nytrue, const int nz, const int nztrue,
              const int ne, const int nJ, const int nw, const int ndwdx,
              const int ndwdp, const int nnz, const int ubw, const int lbw,
              const AMICI_o2mode o2mode, const std::vector<realtype> p,
              const std::vector<realtype> k, const std::vector<int> plist,
              const std::vector<realtype> idlist, const std::vector<int> z2event)
        : nx(nx), nxtrue(nxtrue), ny(ny), nytrue(nytrue),
        nz(nz), nztrue(nztrue), ne(ne), nw(nw), ndwdx(ndwdx), ndwdp(ndwdp),
        nnz(nnz), nJ(nJ), ubw(ubw), lbw(lbw), o2mode(o2mode),
        z2event(z2event),
        idlist(idlist),
        sigmay(ny, 0.0),
        dsigmaydp(ny*plist.size(), 0.0),
        sigmaz(nz, 0.0),
        dsigmazdp(nz*plist.size(), 0.0),
        dJydp(nJ*plist.size(), 0.0),
        dJzdp(nJ*plist.size(), 0.0),
        deltax(nx, 0.0),
        deltasx(nx*plist.size(), 0.0),
        deltaxB(nx, 0.0),
        deltaqB(nJ*plist.size(), 0.0),
        dxdotdp(nx*plist.size(), 0.0),
        x(nx, 0.0),
        sx(plist.size(), std::vector<realtype>(nx, 0.0)),
        y(ny, 0.0),
        my(nytrue, 0.0),
        z(nz, 0.0),
        mz(nztrue, 0.0),
        rz(nz, 0.0),
        dJydy(nJ*nytrue*ny, 0.0),
        dJydsigma(nJ*nytrue*ny, 0.0),
        dJzdz(nJ*nztrue*nz, 0.0),
        dJzdsigma(nJ*nztrue*nz, 0.0),
        dJrzdz(nJ*nztrue*nz, 0.0),
        dJrzdsigma(nJ*nztrue*nz, 0.0),
        dzdx(nz*nx, 0.0),
        dzdp(nz*plist.size(), 0.0),
        drzdx(nz*nx, 0.0),
        drzdp(nz*plist.size(), 0.0),
        dydp(ny*plist.size(), 0.0),
        dydx(ny*nx,0.0),
        w(nw, 0.0),
        dwdx(ndwdx, 0.0),
        dwdp(ndwdp, 0.0),
        M(nx*nx, 0.0),
        stau(plist.size(), 0.0),
        h(ne,0.0),
        p(p),
        k(k),
        plist(plist)
        {
            J = SparseNewMat(nx, nx, nnz, CSC_MAT);
        }
        
        /**
         * @brief Returns a UserData instance with preset model dimensions
         * @return The UserData instance
         */
        UserData getUserData() const;
        
        /**
         * @brief Returns a new UserData instance with preset model dimensions
         * @return The UserData instance
         */
        UserData *getNewUserData() const;
        
        /** Retrieves the solver object
         * @return The Solver instance
         */
        virtual std::unique_ptr<Solver> getSolver() = 0;
        
        /** Root function
         * @param t time
         * @param x state
         * @param dx time derivative of state (DAE only)
         * @param root array to which values of the root function will be written
         */
        virtual void froot(realtype t, AmiVector *x, AmiVector *dx, realtype *root) = 0;
        
        /** Residual function
         * @param t time
         * @param x state
         * @param dx time derivative of state (DAE only)
         * @param xdot array to which values of the residual function will be written
         */
        virtual void fxdot(realtype t, AmiVector *x, AmiVector *dx, AmiVector *xdot) = 0;
        
        /** Dense Jacobian function
         * @param t time
         * @param cj scaling factor (inverse of timestep, DAE only)
         * @param x state
         * @param dx time derivative of state (DAE only)
         * @param xdot values of residual function (unused)
         * @param J dense matrix to which values of the jacobian will be written
         */
        virtual void fJ(realtype t, realtype cj, AmiVector *x, AmiVector *dx,
                              AmiVector *xdot, DlsMat J) = 0;
        
        /** Sparse Jacobian function
         * @param t time
         * @param cj scaling factor (inverse of timestep, DAE only)
         * @param x state
         * @param dx time derivative of state (DAE only)
         * @param xdot values of residual function (unused)
         * @param J sparse matrix to which values of the Jacobian will be written
         */
        virtual void fJSparse(realtype t, realtype cj, AmiVector *x, AmiVector *dx,
                            AmiVector *xdot, SlsMat J) = 0;
        
        /** Diagonal Jacobian function
         * @param t time
         * @param Jdiag array to which the diagonal of the Jacobian will be written
         * @param cj scaling factor (inverse of timestep, DAE only)
         * @param x state
         * @param dx time derivative of state (DAE only)
         * @return flag indicating successful evaluation
         */
        virtual void fJDiag(realtype t, AmiVector *Jdiag, realtype cj, AmiVector *x,
                                AmiVector *dx) = 0;
        
        /** parameter derivative of residual function
         * @param t time
         * @param x state
         * @param dx time derivative of state (DAE only)
         * @return flag indicating successful evaluation
         */
        virtual void fdxdotdp(realtype t, AmiVector *x, AmiVector *dx) = 0;
        
        /** Jacobian multiply function
         * @param t time
         * @param x state
         * @param dx time derivative of state (DAE only)
         * @param xdot values of residual function (unused)
         * @param v multiplication vector (unused)
         * @param nJv array to which result of multiplication will be written
         * @param cj scaling factor (inverse of timestep, DAE only)
         */
        virtual void fJv(realtype t, AmiVector *x, AmiVector *dx, AmiVector *xdot,
                             AmiVector *v, AmiVector *nJv, realtype cj) = 0;
        
        void fx0(AmiVector *x, const UserData *udata);
        
        /** Initial value for time derivative of states (only necessary for DAEs)
         * @param x0 Vector with the initial states
         * @param dx0 Vector to which the initial derivative states will be
         * written (only DAE)
         **/
        virtual void fdx0(AmiVector *x0, AmiVector *dx0) {};

        void fsx0(AmiVectorArray *sx, const AmiVector *x, const UserData *udata);
        
        /** Sensitivity of derivative initial states sensitivities sdx0 (only
         *necessary for DAEs)
         **/
        virtual void fsdx0() {};
        
        void fstau(const realtype t, const int ie, const AmiVector *x, const AmiVectorArray *sx);
        
        void fy(int it, ReturnData *rdata);
        
        void fdydp(const int it, ReturnData *rdata);
        
        void fdydx(const int it, ReturnData *rdata);
        
        void fz(const int nroots, const int ie, const realtype t, const AmiVector *x, ReturnData *rdata);
        
        void fsz(const int nroots, const int ie, const realtype t, const AmiVector *x, const AmiVectorArray *sx, ReturnData *rdata);
        
        void frz(const int nroots, const int ie, const realtype t, const AmiVector *x, ReturnData *rdata);
        
        void fsrz(const int nroots, const int ie, const realtype t, const AmiVector *x, const AmiVectorArray *sx, ReturnData *rdata);
        
        void fdzdp(const realtype t, const int ie, const AmiVector *x);
        
        void fdzdx(const realtype t, const int ie, const AmiVector *x);
        
        void fdrzdp(const realtype t, const int ie, const AmiVector *x);
        
        void fdrzdx(const realtype t, const int ie, const AmiVector *x);
        
        void fdeltax(const int ie, const realtype t, const AmiVector *x,
                             const AmiVector *xdot, const AmiVector *xdot_old);
        
        void fdeltasx(const int ie, const realtype t, const AmiVector *x, const AmiVectorArray *sx,
                              const AmiVector *xdot, const AmiVector *xdot_old);
        
        void fdeltaxB(const int ie, const realtype t, const AmiVector *x, const AmiVector *xB,
                              const AmiVector *xdot, const AmiVector *xdot_old);
        
        void fdeltaqB(const int ie, const realtype t, const AmiVector *x, const AmiVector *xB,
                              const AmiVector *xdot, const AmiVector *xdot_old);
        
        void fsigma_y(const int it, const ExpData *edata, ReturnData *rdata);
        
        void fdsigma_ydp(const int it, const ReturnData *rdata);
        
        void fsigma_z(const realtype t, const int ie, const int *nroots,
                      const ExpData *edata, ReturnData *rdata);
        
        void fdsigma_zdp(const realtype t);
        
        void fJy(const int it, ReturnData *rdata, const ExpData *edata);
        
        void fJz(const int nroots, ReturnData *rdata, const ExpData *edata);
        
        void fJrz(const int nroots, ReturnData *rdata, const ExpData *edata);
        
        void fdJydy(const int it, const ReturnData *rdata, const ExpData *edata);
        
        void fdJydsigma(const int it, const ReturnData *rdata, const ExpData *edata);
        
        void fdJzdz(const int nroots, const ReturnData *rdata, const ExpData *edata);
        
        void fdJzdsigma(const int nroots, const ReturnData *rdata, const ExpData *edata);
        
        void fdJrzdz(const int nroots, const ReturnData *rdata, const ExpData *edata);
        
        void fdJrzdsigma(const int nroots,const ReturnData *rdata, const ExpData *edata);
        
        /** default destructor */
        virtual ~Model() {
            if(J)
                SparseDestroyMat(J);
        };
        
        // Generic implementations
        void fsy(const int it, ReturnData *rdata);
        
        void fsz_tf(const int ie, ReturnData *rdata);
        
        void fsJy(const int it, const std::vector<realtype> dJydx, ReturnData *rdata);
        
        void fdJydp(const int it, const ExpData *edata,
                   const ReturnData *rdata);
        
        void fdJydx(std::vector<realtype> *dJydx, const int it, const ExpData *edata, const ReturnData *rdata);
        
        void fsJz(const int nroots, const std::vector<realtype> dJzdx, AmiVectorArray *sx, const ReturnData *rdata);
        
        void fdJzdp(const int nroots, realtype t, const ExpData *edata, const ReturnData *rdata);
        
        void fdJzdx(std::vector<realtype> *dJzdx, const int nroots, realtype t, const ExpData *edata, const ReturnData *rdata);
        
        void initialize(AmiVector *x, AmiVector *dx, const UserData *udata);
        
        void initializeStates(AmiVector *x, const UserData *udata);
        
        void initHeaviside(AmiVector *x, AmiVector *dx, const UserData *udata);
        
        /** number of paramaeters wrt to which sensitivities are computed
         * @return length of sensitivity index vector
         */
        const int nplist() const {
            return plist.size();
        };
        /** total number of model parameters
         * @return length of parameter vector
         */
        const int np() const {
            return p.size();
        };
        /** number of constants
         * @return length of constant vector
         */
        const int nk() const {
            return k.size();
        };
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
        const std::vector<int> z2event;
        /** flag array for DAE equations */
        const std::vector<realtype> idlist;
        
        /** data standard deviation */
        std::vector<realtype> sigmay;
        /** parameter derivative of data standard deviation */
        std::vector<realtype> dsigmaydp;
        /** event standard deviation */
        std::vector<realtype> sigmaz;
        /** parameter derivative of event standard deviation */
        std::vector<realtype> dsigmazdp;
        /** parameter derivative of data likelihood */
        std::vector<realtype> dJydp;
        /** parameter derivative of event likelihood */
        std::vector<realtype> dJzdp;
        
        /** change in x */
        std::vector<realtype> deltax;
        /** change in sx */
        std::vector<realtype> deltasx;
        /** change in xB */
        std::vector<realtype> deltaxB;
        /** change in qB */
        std::vector<realtype> deltaqB;
        
        /** tempory storage of dxdotdp data across functions */
        std::vector<realtype> dxdotdp;
        
        void fw(const realtype t, const N_Vector x);

        void fdwdp(const realtype t, const N_Vector x);
        
        void fdwdx(const realtype t, const N_Vector x);
        
        /**
         * updateHeaviside updates the heaviside variables h on event occurences
         *
         * @param rootsfound provides the direction of the zero-crossing, so adding
         it will give the right update to the heaviside variables (zero if no root
         was found)
         */
        void updateHeaviside(const std::vector<int> rootsfound) {
            for (int ie = 0; ie < ne; ie++) {
                h.at(ie) += rootsfound.at(ie);
            }
        }
        
        /**
         * updateHeavisideB updates the heaviside variables h on event occurences
         in the backward problem
         *
         * @param rootsfound provides the direction of the zero-crossing, so adding
         it will give the right update to the heaviside variables (zero if no root
         was found)
         */
        void updateHeavisideB(const int *rootsfound) {
            for (int ie = 0; ie < ne; ie++) {
                h.at(ie) -= rootsfound[ie];
            }
        }
        
    protected:
        
        /** model specific implementation of fx0
         * @param x0 initial state
         * @param t initial time
         * @param p parameter vector
         * @param k constant vector
         **/
        virtual void fx0(realtype *x0, const realtype t, const realtype *p, const realtype *k) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fsx0
         * @param sx0 initial state sensitivities
         * @param t initial time
         * @param x0 initial state
         * @param p parameter vector
         * @param k constant vector
         * @param ip sensitivity index
         **/
        virtual void fsx0(realtype *sx0, const realtype t,const realtype *x0, const realtype *p, const realtype *k, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fstau
         * @param stau total derivative of event timepoint
         * @param t current time
         * @param x current state
         * @param p parameter vector
         * @param k constant vector
         * @param h heavyside vector
         * @param sx current state sensitivity
         * @param ip sensitivity index
         * @param ie event index
         **/
        virtual void fstau(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip, const int ie) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fy
         * @param y model output at current timepoint
         * @param t current time
         * @param x current state
         * @param p parameter vector
         * @param k constant vector
         * @param h heavyside vector
         **/
        virtual void fy(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fdydp
         * @param dydp partial derivative of observables y w.r.t. model parameters p
         * @param t current time
         * @param x current state
         * @param p parameter vector
         * @param k constant vector
         * @param h heavyside vector
         * @param ip parameter index w.r.t. which the derivative is requested
         **/
        virtual void fdydp(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fdydx
         * @param dydx partial derivative of observables y w.r.t. model states x
         * @param t current time
         * @param x current state
         * @param p parameter vector
         * @param k constant vector
         * @param h heavyside vector
         **/
        virtual void fdydx(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fz
         * @param z value of event output
         * @param ie event index
         * @param t current time
         * @param x current state
         * @param p parameter vector
         * @param k constant vector
         * @param h heavyside vector
         **/
        virtual void fz(realtype *z, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fsz
         * @param sz Sensitivity of rz, total derivative
         * @param ie event index
         * @param t current time
         * @param x current state
         * @param p parameter vector
         * @param k constant vector
         * @param h heavyside vector
         * @param sx current state sensitivity
         * @param ip sensitivity index
         **/
        virtual void fsz(realtype *sz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of frz
         * @param rz value of root function at current timepoint (non-output events not included)
         * @param ie event index
         * @param t current time
         * @param x current state
         * @param p parameter vector
         * @param k constant vector
         * @param h heavyside vector
         **/
        virtual void frz(realtype *rz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fsrz
         * @param srz Sensitivity of rz, total derivative
         * @param ie event index
         * @param t current time
         * @param x current state
         * @param p parameter vector
         * @param k constant vector
         * @param sx current state sensitivity
         * @param h heavyside vector
         * @param ip sensitivity index
         **/
        virtual void fsrz(realtype *srz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fdzdp
         * @param dzdp partial derivative of event-resolved output z w.r.t. model parameters p
         * @param ie event index
         * @param t current time
         * @param x current state
         * @param p parameter vector
         * @param k constant vector
         * @param h heavyside vector
         * @param ip parameter index w.r.t. which the derivative is requested
         **/
        virtual void fdzdp(realtype *dzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fdzdx
         * @param dzdx partial derivative of event-resolved output z w.r.t. model states x
         * @param ie event index
         * @param t current time
         * @param x current state
         * @param p parameter vector
         * @param k constant vector
         * @param h heavyside vector
         **/
        virtual void fdzdx(realtype *dzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fdzdp
         * @param drzdp partial derivative of root output rz w.r.t. model parameters p
         * @param ie event index
         * @param t current time
         * @param x current state
         * @param p parameter vector
         * @param k constant vector
         * @param h heavyside vector
         * @param ip parameter index w.r.t. which the derivative is requested
         **/
        virtual void fdrzdp(realtype *drzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fdrzdx
         * @param drzdx partial derivative of root output rz w.r.t. model states x
         * @param ie event index
         * @param t current time
         * @param x current state
         * @param p parameter vector
         * @param k constant vector
         * @param h heavyside vector
         **/
        virtual void fdrzdx(realtype *drzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fdeltax
         * @param deltax state update
         * @param t current time
         * @param x current state
         * @param p parameter vector
         * @param k constant vector
         * @param h heavyside vector
         * @param ie event index
         * @param xdot new model right hand side
         * @param xdot_old previous model right hand side
         **/
        virtual void fdeltax(realtype *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h,
                             const int ie, const realtype *xdot, const realtype *xdot_old) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fdeltasx
         * @param deltasx sensitivity update
         * @param t current time
         * @param x current state
         * @param p parameter vector
         * @param k constant vector
         * @param h heavyside vector
         * @param w repeating elements vector
         * @param ip sensitivity index
         * @param ie event index
         * @param xdot new model right hand side
         * @param xdot_old previous model right hand side
         * @param sx state sensitivity
         * @param stau event-time sensitivity
         **/
        virtual void fdeltasx(realtype *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h,
                              const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx,
                              const realtype *stau) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fdeltaxB
         * @param deltaxB adjoint state update
         * @param t current time
         * @param x current state
         * @param p parameter vector
         * @param k constant vector
         * @param h heavyside vector
         * @param ie event index
         * @param xdot new model right hand side
         * @param xdot_old previous model right hand side
         * @param xB current adjoint state
         **/
        virtual void fdeltaxB(realtype *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h,
                              const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fdeltasx
         * @param deltaqB sensitivity update
         * @param t current time
         * @param x current state
         * @param p parameter vector
         * @param k constant vector
         * @param h heavyside vector
         * @param ip sensitivity index
         * @param ie event index
         * @param xdot new model right hand side
         * @param xdot_old previous model right hand side
         * @param xB adjoint state
         **/
        virtual void fdeltaqB(realtype *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h,
                              const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fsigmay
         * @param sigmay standard deviation of measurements
         * @param t current time
         * @param p parameter vector
         * @param k constant vector
         **/
        virtual void fsigma_y(realtype *sigmay, const realtype t, const realtype *p, const realtype *k) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fsigmay
         * @param dsigmaydp partial derivative of standard deviation of measurements
         * @param t current time
         * @param p parameter vector
         * @param k constant vector
         * @param ip sensitivity index
         **/
        virtual void fdsigma_ydp(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fsigmaz
         * @param sigmaz standard deviation of event measurements
         * @param t current time
         * @param p parameter vector
         * @param k constant vector
         **/
        virtual void fsigma_z(realtype *sigmaz, const realtype t, const realtype *p, const realtype *k) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fsigmaz
         * @param dsigmazdp partial derivative of standard deviation of event measurements
         * @param t current time
         * @param p parameter vector
         * @param k constant vector
         * @param ip sensitivity index
         **/
        virtual void fdsigma_zdp(realtype *dsigmazdp, const realtype t, const realtype *p, const realtype *k, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fJy
         * @param nllh negative log-likelihood for measurements y
         * @param iy output index
         * @param p parameter vector
         * @param k constant vector
         * @param y model output at timepoint
         * @param sigmay measurement standard deviation at timepoint
         * @param my measurements at timepoint
         **/
        virtual void fJy(realtype *nllh,const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        /** model specific implementation of fJz
         * @param nllh negative log-likelihood for event measurements z
         * @param iz event output index
         * @param p parameter vector
         * @param k constant vector
         * @param z model event output at timepoint
         * @param sigmaz event measurement standard deviation at timepoint
         * @param mz event measurements at timepoint
         **/
        virtual void fJz(realtype *nllh, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz, const realtype *mz) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fJrz
         * @param nllh regularization for event measurements z
         * @param iz event output index
         * @param p parameter vector
         * @param k constant vector
         * @param z model event output at timepoint
         * @param sigmaz event measurement standard deviation at timepoint
         **/
        virtual void fJrz(realtype *nllh, const int iz, const realtype *p, const realtype *k, const realtype *z, const realtype *sigmaz) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fdJydy
         * @param dJydy partial derivative of time-resolved measurement negative log-likelihood Jy
         * @param iy output index
         * @param p parameter vector
         * @param k constant vector
         * @param y model output at timepoint
         * @param sigmay measurement standard deviation at timepoint
         * @param my measurement at timepoint
         **/
        virtual void fdJydy(realtype *dJydy, const int iy, const realtype *p, const realtype *k,
                            const realtype *y, const realtype *sigmay, const realtype *my) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fdJydsigma
         * @param dJydsigma Sensitivity of time-resolved measurement
         * negative log-likelihood Jy w.r.t. standard deviation sigmay
         * @param iy output index
         * @param p parameter vector
         * @param k constant vector
         * @param y model output at timepoint
         * @param sigmay measurement standard deviation at timepoint
         * @param my measurement at timepoint
         **/
        virtual void fdJydsigma(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k,
                                const realtype *y, const realtype *sigmay, const realtype *my) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fdJzdz
         * @param dJzdz partial derivative of event measurement negative log-likelihood Jz
         * @param iz event output index
         * @param p parameter vector
         * @param k constant vector
         * @param z model event output at timepoint
         * @param sigmaz event measurement standard deviation at timepoint
         * @param mz event measurement at timepoint
         **/
        virtual void fdJzdz(realtype *dJzdz, const int iz, const realtype *p, const realtype *k,
                            const realtype *z, const realtype *sigmaz, const realtype *mz) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fdJzdsigma
         * @param dJzdsigma Sensitivity of event measurement
         * negative log-likelihood Jz w.r.t. standard deviation sigmaz
         * @param iz event output index
         * @param p parameter vector
         * @param k constant vector
         * @param z model event output at timepoint
         * @param sigmaz event measurement standard deviation at timepoint
         * @param mz event measurement at timepoint
         **/
        virtual void fdJzdsigma(realtype *dJzdsigma, const int iz, const realtype *p, const realtype *k,
                                const realtype *z, const realtype *sigmaz, const realtype *mz) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fdJrzdz
         * @param dJrzdz partial derivative of event penalization Jrz
         * @param iz event output index
         * @param p parameter vector
         * @param k constant vector
         * @param rz model root output at timepoint
         * @param sigmaz event measurement standard deviation at timepoint
         **/
        virtual void fdJrzdz(realtype *dJrzdz, const int iz, const realtype *p, const realtype *k,
                             const realtype *rz, const realtype *sigmaz) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fdJrzdsigma
         * @param dJrzdsigma Sensitivity of event penalization Jrz w.r.t.
         * standard deviation sigmaz
         * @param iz event output index
         * @param p parameter vector
         * @param k constant vector
         * @param rz model root output at timepoint
         * @param sigmaz event measurement standard deviation at timepoint
         **/
        virtual void fdJrzdsigma(realtype *dJrzdsigma, const int iz, const realtype *p, const realtype *k,
                                 const realtype *rz, const realtype *sigmaz) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fw
         * @param w Recurring terms in xdot
         * @param t timepoint
         * @param x Vector with the states
         * @param p parameter vector
         * @param k constants vector
         * @param h heavyside vector
         */
        virtual void fw(realtype *w, const realtype t, const realtype *x, const realtype *p,
                        const realtype *k, const realtype *h) {};
        
        /** model specific implementation of dwdp
         * @param dwdp Recurring terms in xdot, parameter derivative
         * @param t timepoint
         * @param x Vector with the states
         * @param p parameter vector
         * @param k constants vector
         * @param h heavyside vector
         * @param w vector with helper variables
         */
        virtual void fdwdp(realtype *dwdp, const realtype t, const realtype *x, const realtype *p,
                           const realtype *k, const realtype *h, const realtype *w) {};
        
        /** model specific implementation of dwdx
         * @param dwdx Recurring terms in xdot, state derivative
         * @param t timepoint
         * @param x Vector with the states
         * @param p parameter vector
         * @param k constants vector
         * @param h heavyside vector
         * @param w vector with helper variables
         */
        virtual void fdwdx(realtype *dwdx, const realtype t, const realtype *x, const realtype *p,
                           const realtype *k, const realtype *h, const realtype *w) {};
        
        void getmy(const int it, const ExpData *edata);
        
        void gety(const int it, const ReturnData *rdata);
        
        void getx(const int it, const ReturnData *rdata);
        
        void getsx(const int it, const ReturnData *rdata);
        
        const realtype gett(const int it, const ReturnData *rdata) const;
        
        void getmz(const int nroots, const ExpData *edata);
        
        void getz(const int nroots, const ReturnData *rdata);
        
        void getrz(const int nroots, const ReturnData *rdata);
        
        /** Jacobian */
        SlsMat J = nullptr;
        
        /** current state */
        std::vector<realtype> x;
        /** current state */
        std::vector<std::vector<realtype>> sx;
        /** current observable */
        std::vector<realtype> y;
        /** current observable measurement */
        std::vector<realtype> my;
        /** current event output */
        std::vector<realtype> z;
        /** current event measurement */
        std::vector<realtype> mz;
        /** current root output */
        std::vector<realtype> rz;
        /** observable derivative of data likelihood */
        std::vector<realtype> dJydy;
        /** observable sigma derivative of data likelihood */
        std::vector<realtype> dJydsigma;
        /** event ouput derivative of event likelihood */
        std::vector<realtype> dJzdz;
        /** event sigma derivative of event likelihood */
        std::vector<realtype> dJzdsigma;
        /** event ouput derivative of event likelihood at final timepoint */
        std::vector<realtype> dJrzdz;
        /** event sigma derivative of event likelihood at final timepoint */
        std::vector<realtype> dJrzdsigma;
        /** state derivative of event output */
        std::vector<realtype> dzdx;
        /** parameter derivative of event output */
        std::vector<realtype> dzdp;
        /** state derivative of event timepoint */
        std::vector<realtype> drzdx;
        /** parameter derivative of event timepoint */
        std::vector<realtype> drzdp;
        /** parameter derivative of observable */
        std::vector<realtype> dydp;
        /** state derivative of observable */
        std::vector<realtype> dydx;
        /** tempory storage of w data across functions */
        std::vector<realtype> w;
        /** tempory storage of dwdx data across functions */
        std::vector<realtype> dwdx;
        /** tempory storage of dwdp data across functions */
        std::vector<realtype> dwdp;
        /** tempory storage of M data across functions */
        std::vector<realtype> M;
        /** tempory storage of stau data across functions */
        std::vector<realtype> stau;
        
        /** flag indicating whether a certain heaviside function should be active or
         not */
        std::vector<realtype> h;
        /** parameters */
        const std::vector<realtype> p;
        /** constants */
        const std::vector<realtype> k;
        /** indexes of parameters wrt to which sensitivities are computed */
        const std::vector<int> plist;
    };
    
} // namespace amici

#endif // fH
