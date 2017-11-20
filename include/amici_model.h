#ifndef AMICI_MODEL_H
#define AMICI_MODEL_H

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
        : nplist(0), np(0), nk(0), nx(0), nxtrue(0), ny(0), nytrue(0), nz(0), nztrue(0),
        ne(0), nw(0), ndwdx(0), ndwdp(0), nnz(0), nJ(0), ubw(0), lbw(0),
        o2mode(AMICI_O2MODE_NONE) {}
        
        /** constructor with model dimensions
         * @param np number of parameters
         * @param nx number of state variables
         * @param nxtrue number of state variables of the non-augmented model
         * @param nk number of constants
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
         */
        Model(const int np, const int nx, const int nxtrue, const int nk,
              const int ny, const int nytrue, const int nz, const int nztrue,
              const int ne, const int nJ, const int nw, const int ndwdx,
              const int ndwdp, const int nnz, const int ubw, const int lbw,
              const AMICI_o2mode o2mode, const std::vector<realtype> p,
              const std::vector<realtype> k, const std::vector<int> plist,
              const std::vector<realtype> idlist, const std::vector<int> z2event)
        : nplist(plist.size()), np(np), nk(nk), nx(nx), nxtrue(nxtrue), ny(ny), nytrue(nytrue),
        nz(nz), nztrue(nztrue), ne(ne), nw(nw), ndwdx(ndwdx), ndwdp(ndwdp),
        nnz(nnz), nJ(nJ), ubw(ubw), lbw(lbw), o2mode(o2mode),
        z2event(z2event),
        idlist(idlist),
        sigmay(ny, 0.0),
        dsigmaydp(ny*nplist, 0.0),
        sigmaz(nz, 0.0),
        dsigmazdp(nz*nplist, 0.0),
        dJydp(nJ*nplist, 0.0),
        dJzdp(nJ*nplist, 0.0),
        deltax(nx, 0.0),
        deltasx(nx*nplist, 0.0),
        deltaxB(nx, 0.0),
        deltaqB(nJ*nplist, 0.0),
        dxdotdp(nx*nplist, 0.0),
        dJydyTmp(nJ * ny, 0.0),
        dJydxTmp(nJ * nx, 0.0),
        dJydsigmaTmp(nJ * ny, 0.0),
        dJzdzTmp(nJ * nz, 0.0),
        dJzdxTmp(nJ * nx, 0.0),
        dJzdsigmaTmp(nJ * nz, 0.0),
        dJrzdsigmaTmp(nJ * nz, 0.0),
        x(nx, 0.0),
        sx(np, std::vector<double>(nx, 0.0)),
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
        dzdp(nz*nplist, 0.0),
        drzdx(nz*nx, 0.0),
        drzdp(nz*nplist, 0.0),
        dydp(ny*nplist, 0.0),
        dydx(ny*nx,0.0),
        w(nw, 0.0),
        dwdx(ndwdx, 0.0),
        dwdp(ndwdp, 0.0),
        M(nx*nx, 0.0),
        stau(nplist, 0.0),
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
         * @return Solver solver object @type Solver
         */
        virtual Solver *getSolver() = 0;
        
        virtual void froot(realtype t, AmiVector *x, AmiVector *dx, realtype *root) = 0;
        
        virtual void fxdot(realtype t, AmiVector *x, AmiVector *dx, AmiVector *xdot) = 0;
        
        virtual void fJ(realtype t, realtype cj, AmiVector *x, AmiVector *dx,
                              AmiVector *xdot, DlsMat J) = 0;
        
        virtual void fJSparse(realtype t, realtype cj, AmiVector *x, AmiVector *dx,
                            AmiVector *xdot, SlsMat J) = 0;
        
        virtual int fJDiag(realtype t, AmiVector *Jdiag, realtype cj, AmiVector *x,
                                AmiVector *dx) = 0;
        
        virtual int fdxdotdp(realtype t, AmiVector *x, AmiVector *dx) = 0;
        
        virtual void fJv(realtype t, AmiVector *x, AmiVector *dx, AmiVector *xdot,
                             AmiVector *v, AmiVector *nJv, realtype cj) = 0;
        

        void fx0(AmiVector *x, const UserData *udata);
        
        /** model specific implementation of fx0
         * param sx0 initial state sensitivities
         * param t initial time
         * param x0 initial state
         * param p parameter vector
         * param k constant vector
         * param ip parameter index
         **/
        virtual void model_x0(realtype *x0, const realtype t, const realtype *p, const realtype *k) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** Initial value for time derivative of states (only necessary for DAEs)
         * @param x0 Vector with the initial states @type N_Vector
         * @param dx0 Vector to which the initial derivative states will be
         *written (only DAE) @type N_Vector
         * @param user_data object with model specifications @type TempData
         **/
        virtual void fdx0(AmiVector *x0, AmiVector *dx0) {};
        
        void fsx0(AmiVectorArray *sx, const AmiVector *x, const UserData *udata);
        
        /** model specific implementation of fsx0
         * param sx0 initial state sensitivities
         * param t initial time
         * param x0 initial state
         * param p parameter vector
         * param k constant vector
         * param ip sensitivity index
         **/
        virtual void model_sx0(realtype *sx0, const realtype t,const realtype *x0, const realtype *p, const realtype *k, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        
        /** Sensitivity of derivative initial states sensitivities sdx0 (only
         *necessary for DAEs)
         * @param udata object with user input
         **/
        virtual void fsdx0() {};
        
        void fstau(const realtype t, const int ie, const AmiVector *x, const AmiVectorArray *sx);
        
        /** model specific implementation of fstau
         * param stau total derivative of event timepoint
         * param t current time
         * param x current state
         * param p parameter vector
         * param k constant vector
         * param sx current state sensitivity
         * param ip sensitivity index
         * param ie event index
         **/
        virtual void model_stau(double *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip, const int ie) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fy(int it, ReturnData *rdata);
        
        /** model specific implementation of fy
         * param y model output at current timepoint
         * param t current time
         * param x current state
         * param p parameter vector
         * param k constant vector
         **/
        virtual void model_y(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }

        void fdydp(const int it, ReturnData *rdata);
        
        /** model specific implementation of fdydp
         * param dydp partial derivative of observables y w.r.t. model parameters p
         * param t current time
         * param x current state
         * param p parameter vector
         * param k constant vector
         * param ip parameter index w.r.t. which the derivative is requested
         **/
        virtual void model_dydp(double *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fdydx(const int it, ReturnData *rdata);
        
        /** model specific implementation of fdydx
         * param dydx partial derivative of observables y w.r.t. model states x
         * param t current time
         * param x current state
         * param p parameter vector
         * param k constant vector
         **/
        virtual void model_dydx(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fz(const int nroots, const int ie, const realtype t, const AmiVector *x, ReturnData *rdata);
        
        /** model specific implementation of fz
         * param z value of event output
         * param t current time
         * param x current state
         * param p parameter vector
         * param k constant vector
         **/
        virtual void model_z(double *z, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fsz(const int nroots, const int ie, const realtype t, const AmiVector *x, const AmiVectorArray *sx, ReturnData *rdata);
        
        /** model specific implementation of fsz
         * param sz Sensitivity of rz, total derivative
         * param t current time
         * param x current state
         * param p parameter vector
         * param k constant vector
         * param sx current state sensitivity
         * param ip sensitivity index
         **/
        virtual void model_sz(double *sz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void frz(const int nroots, const int ie, const realtype t, const AmiVector *x, ReturnData *rdata);
        
        /** model specific implementation of frz
         * param rz value of root function at current timepoint (non-output events not included)
         * param t current time
         * param x current state
         * param p parameter vector
         * param k constant vector
         **/
        virtual void model_rz(double *rz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fsrz(const int nroots, const int ie, const realtype t, const AmiVector *x, const AmiVectorArray *sx, ReturnData *rdata);
        
        /** model specific implementation of fsrz
         * param srz Sensitivity of rz, total derivative
         * param t current time
         * param x current state
         * param p parameter vector
         * param k constant vector
         * param sx current state sensitivity
         * param ip sensitivity index
         **/
        virtual void model_srz(double *srz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fdzdp(const realtype t, const int ie, const AmiVector *x);
        
        /** model specific implementation of fdzdp
         * param dzdp partial derivative of event-resolved output z w.r.t. model parameters p
         * param t current time
         * param x current state
         * param p parameter vector
         * param k constant vector
         * param ip parameter index w.r.t. which the derivative is requested
         **/
        virtual void model_dzdp(double *dzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fdzdx(const realtype t, const int ie, const AmiVector *x);
        
        /** model specific implementation of fdzdx
         * param dzdx partial derivative of event-resolved output z w.r.t. model states x
         * param t current time
         * param x current state
         * param p parameter vector
         * param k constant vector
         **/
        virtual void model_dzdx(double *dzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fdrzdp(const realtype t, const int ie, const AmiVector *x);
        
        /** model specific implementation of fdzdp
         * param drzdp partial derivative of root output rz w.r.t. model parameters p
         * param t current time
         * param x current state
         * param p parameter vector
         * param k constant vector
         * param ip parameter index w.r.t. which the derivative is requested
         **/
        virtual void model_drzdp(double *drzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fdrzdx(const realtype t, const int ie, const AmiVector *x);
        
        /** model specific implementation of fdrzdx
         * param drzdx partial derivative of root output rz w.r.t. model states x
         * param t current time
         * param x current state
         * param p parameter vector
         * param k constant vector
         **/
        virtual void model_drzdx(double *drzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fdeltax(const int ie, const realtype t, const AmiVector *x,
                             const AmiVector *xdot, const AmiVector *xdot_old);
        
        /** model specific implementation of fdeltax
         * param deltax state update
         * param t current time
         * param x current state
         * param p parameter vector
         * param k constant vector
         * param ie event index
         * param xdot new model right hand side
         * param xdot_old previous model right hand side
         **/
        virtual void model_deltax(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h,
                                 const int ie, const realtype *xdot, const realtype *xdot_old) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fdeltasx(const int ie, const realtype t, const AmiVector *x, const AmiVectorArray *sx,
                              const AmiVector *xdot, const AmiVector *xdot_old);
        
        /** model specific implementation of fdeltasx
         * param deltasx sensitivity update
         * param t current time
         * param x current state
         * param p parameter vector
         * param k constant vector
         * param ip sensitivity index
         * param ie event index
         * param xdot new model right hand side
         * param xdot_old previous model right hand side
         * param sx state sensitivity
         * param stau event-time sensitivity
         **/
        virtual void model_deltasx(double *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h,
                                   const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx,
                                   const realtype *stau) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fdeltaxB(const int ie, const realtype t, const AmiVector *x, const AmiVector *xB,
                              const AmiVector *xdot, const AmiVector *xdot_old);
        
        /** model specific implementation of fdeltaxB
         * param deltaxB adjoint state update
         * param t current time
         * param x current state
         * param p parameter vector
         * param k constant vector
         * param ie event index
         * param xdot new model right hand side
         * param xdot_old previous model right hand side
         * param xB current adjoint state
         **/
        virtual void model_deltaxB(double *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h,
                                  const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fdeltaqB(const int ie, const realtype t, const AmiVector *x, const AmiVector *xB,
                              const AmiVector *xdot, const AmiVector *xdot_old);
        
        /** model specific implementation of fdeltasx
         * param deltasx sensitivity update
         * param t current time
         * param x current state
         * param p parameter vector
         * param k constant vector
         * param ip sensitivity index
         * param ie event index
         * param xdot new model right hand side
         * param xdot_old previous model right hand side
         * param xB adjoint state
         * param qBdot right hand side of quadradature
         **/
        virtual void model_deltaqB(double *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h,
                                   const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }

        void fsigma_y(const int it, const ExpData *edata, ReturnData *rdata);
        
        /** model specific implementation of fsigmay
         * param sigmay standard deviation of measurements
         * param t current time
         * param p parameter vector
         * param k constant vector
         **/
        virtual void model_sigma_y(double *sigmay, const realtype t, const realtype *p, const realtype *k) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fdsigma_ydp(const int it, const ReturnData *rdata);
        
        /** model specific implementation of fsigmay
         * param dsigmaydp partial derivative of standard deviation of measurements
         * param t current time
         * param p parameter vector
         * param k constant vector
         **/
        virtual void model_dsigma_ydp(double *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fsigma_z(const realtype t, const int ie, const int *nroots,
                      const ExpData *edata, ReturnData *rdata);
        
        /** model specific implementation of fsigmaz
         * param sigmaz standard deviation of event measurements
         * param t current time
         * param p parameter vector
         * param k constant vector
         **/
        virtual void model_sigma_z(double *sigmaz, const realtype t, const realtype *p, const realtype *k) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fdsigma_zdp(const realtype t);
        
        /** model specific implementation of fsigmaz
         * param dsigmazdp partial derivative of standard deviation of event measurements
         * param t current time
         * param p parameter vector
         * param k constant vector
         **/
        virtual void model_dsigma_zdp(double *dsigmazdp, const realtype t, const realtype *p, const realtype *k, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fJy(const int it, ReturnData *rdata, const ExpData *edata);
        
        /** model specific implementation of fJy
         * param nllh negative log-likelihood for measurements y
         * param p parameter vector
         * param k constant vector
         * param y model output at timepoint
         * param sigmay measurement standard deviation at timepoint
         * param my measurements at timepoint
         **/
        virtual void model_Jy(double *nllh,const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fJz(const int nroots, ReturnData *rdata, const ExpData *edata);
        
        /** model specific implementation of fJz
         * param nllh negative log-likelihood for event measurements z
         * param p parameter vector
         * param k constant vector
         * param z model event output at timepoint
         * param sigmaz event measurement standard deviation at timepoint
         * param mz event measurements at timepoint
         **/
        virtual void model_Jz(double *nllh, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fJrz(const int nroots, ReturnData *rdata, const ExpData *edata);
        
        /** model specific implementation of fJrz
         * param nllh regularization for event measurements z
         * param p parameter vector
         * param k constant vector
         * param z model event output at timepoint
         * param sigmaz event measurement standard deviation at timepoint
         **/
        virtual void model_Jrz(double *nllh, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fdJydy(const int it, const ReturnData *rdata, const ExpData *edata);
        
        /** model specific implementation of fdJydy
         * param dJydy partial derivative of time-resolved measurement negative log-likelihood Jy
         * param p parameter vector
         * param k constant vector
         * param y model output at timepoint
         * param sigmay measurement standard deviation at timepoint
         * param my measurement at timepoint
         **/
        virtual void model_dJydy(double *dJydy, const int iy, const realtype *p, const realtype *k,
                                 const double *y, const double *sigmay, const double *my) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fdJydsigma(const int it, const ReturnData *rdata, const ExpData *edata);
        
        /** model specific implementation of fdJydsigma
         * param dJydsigma Sensitivity of time-resolved measurement
         * negative log-likelihood Jy w.r.t. standard deviation sigmay
         * param p parameter vector
         * param k constant vector
         * param y model output at timepoint
         * param sigmaz measurement standard deviation at timepoint
         * param my measurement at timepoint
         **/
        virtual void model_dJydsigma(double *dJydsigma, const int iy, const realtype *p, const realtype *k,
                                 const double *y, const double *sigmay, const double *my) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fdJzdz(const int nroots, const ReturnData *rdata, const ExpData *edata);
        
        /** model specific implementation of fdJzdz
         * param dJzdz partial derivative of event measurement negative log-likelihood Jz
         * param p parameter vector
         * param k constant vector
         * param z model event output at timepoint
         * param sigmaz event measurement standard deviation at timepoint
         * param mz event measurement at timepoint
         **/
        virtual void model_dJzdz(double *dJzdz, const int iz, const realtype *p, const realtype *k,
                                 const double *z, const double *sigmaz, const double *mz) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fdJzdsigma(const int nroots, const ReturnData *rdata, const ExpData *edata);
        
        /** model specific implementation of fdJzdsigma
         * param dJzdsigma Sensitivity of event measurement
         * negative log-likelihood Jz w.r.t. standard deviation sigmaz
         * param p parameter vector
         * param k constant vector
         * param z model event output at timepoint
         * param sigmaz event measurement standard deviation at timepoint
         * param mz event measurement at timepoint
         **/
        virtual void model_dJzdsigma(double *dJzdsigma, const int iz, const realtype *p, const realtype *k,
                                     const double *z, const double *sigmaz, const double *mz) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fdJrzdz(const int nroots, const ReturnData *rdata, const ExpData *edata);
        
        /** model specific implementation of fdJrzdz
         * param dJrzdz partial derivative of event penalization Jrz
         * param p parameter vector
         * param k constant vector
         * param rz model root output at timepoint
         * param sigmaz event measurement standard deviation at timepoint
         **/
        virtual void model_dJrzdz(double *dJrzdz, const int iz, const realtype *p, const realtype *k,
                                 const double *rz, const double *sigmaz) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        void fdJrzdsigma(const int nroots,const ReturnData *rdata, const ExpData *edata);
        
        /** model specific implementation of fdJrzdsigma
         * param dJzdsigma Sensitivity of event penalization Jz w.r.t.
         * standard deviation sigmaz
         * param p parameter vector
         * param k constant vector
         * param rz model root output at timepoint
         * param sigmaz event measurement standard deviation at timepoint
         **/
        virtual void model_dJrzdsigma(double *dJrzdsigma, const int iz, const realtype *p, const realtype *k,
                                     const double *rz, const double *sigmaz) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        virtual ~Model() {
            if(J)
                SparseDestroyMat(J);
        };
        
        // Generic implementations
        void fsy(const int it, ReturnData *rdata);
        
        void fsz_tf(const int ie, ReturnData *rdata);
        
        void fsJy(const int it, const std::vector<double> dJydx, ReturnData *rdata);
        
        void fdJydp(const int it, const ExpData *edata,
                   const ReturnData *rdata);
        
        void fdJydx(std::vector<double> dJydx, const int it, const ExpData *edata, const ReturnData *rdata);
        
        void fsJz(const int nroots, const std::vector<double> dJzdx, AmiVectorArray *sx, const ReturnData *rdata);
        
        void fdJzdp(const int nroots, realtype t, const ExpData *edata, const ReturnData *rdata);
        
        void fdJzdx(std::vector<double> dJzdx, const int nroots, realtype t, const ExpData *edata, const ReturnData *rdata);
        
        void initialize(AmiVector *x, AmiVector *dx, const UserData *udata);
        
        void initializeStates(AmiVector *x, const UserData *udata);
        
        void initHeaviside(AmiVector *x, AmiVector *dx, const UserData *udata);
        
        // Model dimensions 
        const int nplist;
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
        const std::vector<int> z2event;
        /** flag array for DAE equations */
        const std::vector<realtype> idlist;
        
        /** data standard deviation */
        std::vector<double> sigmay;
        /** parameter derivative of data standard deviation */
        std::vector<double> dsigmaydp;
        /** event standard deviation */
        std::vector<double> sigmaz;
        /** parameter derivative of event standard deviation */
        std::vector<double> dsigmazdp;
        /** parameter derivative of data likelihood */
        std::vector<double> dJydp;
        /** parameter derivative of event likelihood */
        std::vector<double> dJzdp;
        
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
        
        void updateHeaviside(const std::vector<int> rootsfound) {
            /**
             * updateHeaviside updates the heaviside variables h on event occurences
             *
             * @param rootsfound provides the direction of the zero-crossing, so adding
             it will give the right update to the heaviside variables (zero if no root
             was found)
             */
            for (int ie = 0; ie < ne; ie++) {
                h.at(ie) += rootsfound.at(ie);
            }
        }
        
        void updateHeavisideB(const int *rootsfound) {
            /**
             * updateHeavisideB updates the heaviside variables h on event occurences
             in the backward problem
             *
             * @param rootsfound provides the direction of the zero-crossing, so adding
             it will give the right update to the heaviside variables (zero if no root
             was found)
             */
            for (int ie = 0; ie < ne; ie++) {
                h.at(ie) -= rootsfound[ie];
            }
        }
        
        void fw(const realtype t, const N_Vector x);
        
        /** model specific implementation of fw
         * @param w Recurring terms in xdot
         * @param t timepoint
         * @param x Vector with the states
         * @param p parameter vector
         * @param k constants vector
         */
        virtual void model_w(realtype *w, const realtype t, const realtype *x, const realtype *p,
                                const realtype *k, const realtype *h) {};
        
        void fdwdp(const realtype t, const N_Vector x);
        
        /** model specific implementation of dwdp
         * @param dwdp Recurring terms in xdot, parameter derivative
         * @param t timepoint
         * @param x Vector with the states
         * @param p parameter vector
         * @param k constants vector
         * @param w vector with helper variables
         */
        virtual void model_dwdp(realtype *dwdp, const realtype t, const realtype *x, const realtype *p,
                                const realtype *k, const realtype *h, const realtype *w) {};
        
        void fdwdx(const realtype t, const N_Vector x);
        
        /** model specific implementation of dwdx
         * @param dwdx Recurring terms in xdot, state derivative
         * @param t timepoint
         * @param x Vector with the states
         * @param p parameter vector
         * @param k constants vector
         * @param w vector with helper variables
         */
        virtual void model_dwdx(realtype *dwdx, const realtype t, const realtype *x, const realtype *p,
                                const realtype *k, const realtype *h, const realtype *w) {};
        
    protected:
        
        void getmy(const int it, const ExpData *edata);
        
        void gety(const int it, const ReturnData *rdata);
        
        void getx(const int it, const ReturnData *rdata);
        
        void getsx(const int it, const ReturnData *rdata);
        
        const realtype gett(const int it, const ReturnData *rdata) const;
        
        void getmz(const int nroots, const ExpData *edata);
        
        void getz(const int nroots, const ReturnData *rdata);
        
        void getrz(const int nroots, const ReturnData *rdata);

        /** storage for dJydy slice */
        std::vector<double> dJydyTmp;
        /** storage for dJydx slice */
        std::vector<double> dJydxTmp;
        /** storage for dJydsigma slice */
        std::vector<double> dJydsigmaTmp;
        /** storage for dJzdz slice */
        std::vector<double> dJzdzTmp;
        /** storage for dJzdx slice */
        std::vector<double> dJzdxTmp;
        /** storage for dJzdsigma slice */
        std::vector<double> dJzdsigmaTmp;
        /** storage for dJrzdsigma slice */
        std::vector<double> dJrzdsigmaTmp;
        
        /** Jacobian */
        SlsMat J = nullptr;
        
        /** current state */
        std::vector<double> x;
        /** current state */
        std::vector<std::vector<double>> sx;
        /** current observable */
        std::vector<double> y;
        /** current observable measurement */
        std::vector<double> my;
        /** current event output */
        std::vector<double> z;
        /** current event measurement */
        std::vector<double> mz;
        /** current root output */
        std::vector<double> rz;
        /** observable derivative of data likelihood */
        std::vector<double> dJydy;
        /** observable sigma derivative of data likelihood */
        std::vector<double> dJydsigma;
        /** event ouput derivative of event likelihood */
        std::vector<double> dJzdz;
        /** event sigma derivative of event likelihood */
        std::vector<double> dJzdsigma;
        /** event ouput derivative of event likelihood at final timepoint */
        std::vector<double> dJrzdz;
        /** event sigma derivative of event likelihood at final timepoint */
        std::vector<double> dJrzdsigma;
        /** state derivative of event output */
        std::vector<double> dzdx;
        /** parameter derivative of event output */
        std::vector<double> dzdp;
        /** state derivative of event timepoint */
        std::vector<double> drzdx;
        /** parameter derivative of event timepoint */
        std::vector<double> drzdp;
        /** parameter derivative of observable */
        std::vector<double> dydp;
        /** state derivative of observable */
        std::vector<double> dydx;
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
        const std::vector<realtype> p;
        const std::vector<realtype> k;
        const std::vector<int> plist;
        
        friend class IDASolver;
        friend class CVodeSolver;
    };
    
} // namespace amici

#endif // MODEL_H
