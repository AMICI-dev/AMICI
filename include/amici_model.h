#ifndef AMICI_MODEL_H
#define AMICI_MODEL_H


#include <include/rdata.h>
#include <include/udata.h>
#include <include/edata.h>
#include <include/amici.h>
#include <include/amici_exception.h>
#include <include/amici_vector.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_direct.h>
#include <sundials/sundials_sparse.h>

#include <vector>

namespace amici {
    
    /**
     * @brief The Model class represents an AMICI ODE model.
     * The model can compute various model related quantities based
     * on symbolically generated code.
     */
    class Model {
    public:
        /** default constructor */
        Model()
        : np(0), nk(0), nx(0), nxtrue(0), ny(0), nytrue(0), nz(0), nztrue(0),
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
              const AMICI_o2mode o2mode)
        : np(np), nk(nk), nx(nx), nxtrue(nxtrue), ny(ny), nytrue(nytrue),
        nz(nz), nztrue(nztrue), ne(ne), nw(nw), ndwdx(ndwdx), ndwdp(ndwdp),
        nnz(nnz), nJ(nJ), ubw(ubw), lbw(lbw), o2mode(o2mode), nplist(np),
        dJydy(nJ*nytrue*ny, 0.0),
        dJydsigma(nJ*nytrue*ny, 0.0),
        dJzdz(nJ*nztrue*nz, 0.0),
        dJzdsigma(nJ*nztrue*nz, 0.0),
        dJrzdz(nJ*nztrue*nz, 0.0),
        dJrzdsigma(nJ*nztrue*nz, 0.0),
        dJydyTmp(nJ * ny, 0.0),
        dJydxTmp(nJ * nx, 0.0),
        dJydsigmaTmp(nJ * ny, 0.0),
        dJzdzTmp(nJ * nz, 0.0),
        dJzdxTmp(nJ * nx, 0.0),
        dJzdsigmaTmp(nJ * nz, 0.0),
        dJrzdsigmaTmp(nJ * nz, 0.0),
        dzdx(nz*nx, 0.0),
        dzdp(nz*nplist, 0.0),
        drzdx(nz*nx, 0.0),
        drzdp(nz*nplist, 0.0),
        dydp(ny*nplist, 0.0),
        dydx(ny*nx,0.0),
        sigmay(ny, 0.0),
        dsigmaydp(ny*nplist, 0.0),
        sigmaz(nz, 0.0),
        dsigmazdp(nz*nplist, 0.0),
        dxdotdp(nx*nplist, 0.0),
        w(nw, 0.0),
        dwdx(ndwdx, 0.0),
        dwdp(ndwdp, 0.0),
        M(nx*nx, 0.0),
        stau(nplist, 0.0),
        deltax(nx, 0.0),
        deltasx(nx*nplist, 0.0),
        deltaxB(nx, 0.0),
        deltaqB(nJ*nplist, 0.0)
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
        virtual Solver *getSolver() { return NULL; }
        
        /** Initial states
         * @param[in] udata object with user input
         **/
        virtual void fx0(AmiVector x, UserData *udata) {
            x.reset();
            model_x0(x.data(),udata->tstart(),udata->p(),udata->k());
        };
        
        /** model specific implementation of fx0
         * param[out] sx0 initial state sensitivities
         * param[in] t initial time
         * param[in] x0 initial state
         * param[in] p parameter vector
         * param[in] k constant vector
         * param[in] ip parameter index
         **/
        virtual void model_x0(realtype *x0, const realtype t, const realtype *p, const realtype *k) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** Initial value for time derivative of states (only necessary for DAEs)
         * @param[in] x0 Vector with the initial states @type N_Vector
         * @param[out] dx0 Vector to which the initial derivative states will be
         *written (only DAE) @type N_Vector
         * @param[in] user_data object with model specifications @type TempData
         **/
        virtual void fdx0(N_Vector x0, N_Vector dx0, void *user_data) = 0;
        
        /** Initial value for initial state sensitivities
         * @param[in] udata object with user input
         **/
        virtual void fsx0(AmiVectorArray sx, const AmiVector x, const UserData *udata) {
            sx.reset();
            for(int ip = 0; ip<udata->nplist(); ip++)
                model_sx0(sx.data(ip),udata->tstart(),x.data(),udata->p(),udata->k(),udata->plist(ip));
        }
        
        /** model specific implementation of fsx0
         * param[out] sx0 initial state sensitivities
         * param[in] t initial time
         * param[in] x0 initial state
         * param[in] p parameter vector
         * param[in] k constant vector
         * param[in] ip sensitivity index
         **/
        virtual void model_sx0(realtype *sx0, const realtype t,const realtype *x0, const realtype *p, const realtype *k, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        
        /** Sensitivity of derivative initial states sensitivities sdx0 (only
         *necessary for DAEs)
         * @param[in] udata object with user input
         **/
        virtual void fsdx0(const UserData *udata) = 0;
        
        /** Sensitivity of event timepoint, total derivative
         * @param[in] ie event index
         * @param[in] udata object with user input
         */
        virtual void fstau(const int ie,const realtype t, const AmiVector x, const AmiVectorArray sx, const UserData *udata) {
            std::fill(stau.begin(),stau.end(),0.0);
            for(int ip = 0; ip < udata->nplist(); ip++){
                model_stau(stau.data(),t,x.data(),udata->p(),udata->k(),sx.data(ip),udata->plist(ip),ie);
            }
        }
        
        /** model specific implementation of fstau
         * param[out] stau total derivative of event timepoint
         * param[in] t current time
         * param[in] x current state
         * param[in] p parameter vector
         * param[in] k constant vector
         * param[in] sx current state sensitivity
         * param[in] ip sensitivity index
         * param[in] ie event index
         **/
        virtual void model_stau(double *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *sx, const int ip, const int ie) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        
        /** Observables / measurements
         * @param[in] it timepoint index
         * @param[out] rdata pointer to return data object
         * @param[in] udata object with user input
         */
        virtual void fy(int it, const realtype t, const AmiVector x, ReturnData *rdata,
                        const UserData *udata) {
            std::vector<double> yreturn(ny,0.0);
            model_y(yreturn.data(),t,x.data(),udata->p(),udata->k());
            for(int iy; iy < ny; iy++)
                rdata->y[it + udata->nt()*iy] = yreturn.at(iy);
        }
        
        /** model specific implementation of fy
         * param[out] y model output at current timepoint
         * param[in] t current time
         * param[in] x current state
         * param[in] p parameter vector
         * param[in] k constant vector
         **/
        virtual void model_y(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** partial derivative of observables y w.r.t. model parameters p
         * @param[in] udata object with user input
         */
        virtual void fdydp(const realtype t, const AmiVector x, const UserData *udata) {
            std::fill(dydp.begin(),dydp.end(),0.0);
            for(int ip = 0; ip < udata->nplist(); ip++){
                model_dydp(dydp.data(),t,x.data(),udata->p(),udata->k(),udata->plist(ip));
            }
        }
        
        /** model specific implementation of fdydp
         * param[out] dydp partial derivative of observables y w.r.t. model parameters p
         * param[in] t current time
         * param[in] x current state
         * param[in] p parameter vector
         * param[in] k constant vector
         * param[in] ip parameter index w.r.t. which the derivative is requested
         **/
        virtual void model_dydp(double *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** partial derivative of observables y w.r.t. state variables x
         const UserData *udata
         */
        virtual void fdydx(const realtype t, const AmiVector x, const UserData *udata) {
            std::fill(dydx.begin(),dydx.end(),0.0);
            model_dydx(dydx.data(),t,x.data(),udata->p(),udata->k());
        }
        
        /** model specific implementation of fdydx
         * param[out] dydx partial derivative of observables y w.r.t. model states x
         * param[in] t current time
         * param[in] x current state
         * param[in] p parameter vector
         * param[in] k constant vector
         **/
        virtual void model_dydx(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** Event-resolved output
         * @param[in] nroots number of events for event index
         * @param[out] rdata pointer to return data object
         * @param[in] udata object with user input
         */
        virtual void fz(const int nroots, const realtype t, const AmiVector x, ReturnData *rdata,
                        const UserData *udata) {
            std::vector<double> zreturn(nz,0.0);
            model_z(zreturn.data(),t,x.data(),udata->p(),udata->k());
            for(int iz; iz < nz; iz++) {
                rdata->z[nroots+udata->nme()*iz] = zreturn.at(iz);
            }
        }
        
        /** model specific implementation of fz
         * param[out] z value of event output
         * param[in] t current time
         * param[in] x current state
         * param[in] p parameter vector
         * param[in] k constant vector
         **/
        virtual void model_z(double *z, const realtype t, const realtype *x, const realtype *p, const realtype *k) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** Sensitivity of z, total derivative
         * @param[in] nroots number of events for event index
         * @param[out] rdata pointer to return data object
         * @param[in] udata object with user input
         */
        virtual void fsz(const int nroots, const realtype t, const AmiVector x, const AmiVectorArray sx, ReturnData *rdata,
                         const UserData *udata) {
            for(int ip; ip < udata->nplist();  ip++ ){
                std::vector<double> szreturn(nz,0.0);
                model_sz(szreturn.data(),t,x.data(),udata->p(),udata->k(),sx.data(ip),udata->plist(ip));
                for(int iz; iz < nz; iz++) {
                    rdata->sz[nroots+udata->nme()*(ip*nz + iz)] = szreturn.at(iz);
                }
            }
        }
        
        /** model specific implementation of fsz
         * param[out] sz Sensitivity of rz, total derivative
         * param[in] t current time
         * param[in] x current state
         * param[in] p parameter vector
         * param[in] k constant vector
         * param[in] sx current state sensitivity
         * param[in] ip sensitivity index
         **/
        virtual void model_sz(double *sz, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *sx, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** Event root function of events (equal to froot but does not include
         * non-output events)
         * @param[in] nroots number of events for event index
         * @param[in] udata object with user input
         * @param[out] rdata pointer to return data object
         */
        virtual void frz(const int nroots, const realtype t, const AmiVector x, ReturnData *rdata,
                         const UserData *udata) {
            std::vector<double> rzreturn(nz,0.0);
            model_rz(rzreturn.data(),t,x.data(),udata->p(),udata->k());
            for(int iz; iz < nz; iz++) {
                rdata->rz[nroots+udata->nme()*iz] = rzreturn.at(iz);
            }
        }
        
        /** model specific implementation of frz
         * param[out] rz value of root function at current timepoint (non-output events not included)
         * param[in] t current time
         * param[in] x current state
         * param[in] p parameter vector
         * param[in] k constant vector
         **/
        virtual void model_rz(double *rz, const realtype t, const realtype *x, const realtype *p, const realtype *k) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** Sensitivity of rz, total derivative
         * @param[in] nroots number of events for event index
         * @param[out] rdata pointer to return data object
         * @param[in] udata object with user input
         */
        virtual void fsrz(const int nroots, const realtype t, const AmiVector x, const AmiVectorArray sx, ReturnData *rdata,
                          const UserData *udata) {
            for(int ip; ip < udata->nplist();  ip++ ){
                std::vector<double> srzreturn(nz,0.0);
                model_srz(srzreturn.data(),t,x.data(),udata->p(),udata->k(),sx.data(ip),udata->plist(ip));
                for(int iz; iz < nz; iz++) {
                    rdata->srz[nroots+udata->nme()*(ip*nz + iz)] = srzreturn.at(iz);
                }
            }
        }
        
        /** model specific implementation of fsrz
         * param[out] srz Sensitivity of rz, total derivative
         * param[in] t current time
         * param[in] x current state
         * param[in] p parameter vector
         * param[in] k constant vector
         * param[in] sx current state sensitivity
         * param[in] ip sensitivity index
         **/
        virtual void model_srz(double *srz, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *sx, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** partial derivative of event-resolved output z w.r.t. to model parameters p
         * @param[in] udata object with user input
         */
        virtual void fdzdp(const realtype t, const AmiVector x, const UserData *udata) {
            std::fill(dzdp.begin(),dzdp.end(),0.0);
            for(int ip = 0; ip < udata->nplist(); ip++){
                model_dzdp(dzdp.data(),t,x.data(),udata->p(),udata->k(),udata->plist(ip));
            }
        }
        
        /** model specific implementation of fdzdp
         * param[out] dzdp partial derivative of event-resolved output z w.r.t. model parameters p
         * param[in] t current time
         * param[in] x current state
         * param[in] p parameter vector
         * param[in] k constant vector
         * param[in] ip parameter index w.r.t. which the derivative is requested
         **/
        virtual void model_dzdp(double *dzdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** partial derivative of event-resolved output z w.r.t. to model states x
         * @param[in] udata object with user input
         */
        virtual void fdzdx(const realtype t, const AmiVector x, const UserData *udata) {
            std::fill(dzdx.begin(),dzdx.end(),0.0);
            model_dzdx(dzdx.data(),t,x.data(),udata->p(),udata->k());
        }
        
        /** model specific implementation of fdzdx
         * param[out] dzdx partial derivative of event-resolved output z w.r.t. model states x
         * param[in] t current time
         * param[in] x current state
         * param[in] p parameter vector
         * param[in] k constant vector
         **/
        virtual void model_dzdx(double *dzdx, const realtype t, const realtype *x, const realtype *p, const realtype *k) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** Sensitivity of event-resolved root output w.r.t. to model parameters p
         * @param[in] udata object with user input
         */
        virtual void fdrzdp(const realtype t, const AmiVector x, const UserData *udata) {
            std::fill(drzdp.begin(),drzdp.end(),0.0);
            for(int ip = 0; ip < udata->nplist(); ip++){
                model_drzdp(drzdp.data(),t,x.data(),udata->p(),udata->k(),udata->plist(ip));
            }
        }
        
        /** model specific implementation of fdzdp
         * param[out] drzdp partial derivative of root output rz w.r.t. model parameters p
         * param[in] t current time
         * param[in] x current state
         * param[in] p parameter vector
         * param[in] k constant vector
         * param[in] ip parameter index w.r.t. which the derivative is requested
         **/
        virtual void model_drzdp(double *drzdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** Sensitivity of event-resolved measurements rz w.r.t. to model states x
         * @param[in] udata object with user input
         */
        virtual void fdrzdx(const realtype t, const AmiVector x, const UserData *udata) {
            std::fill(drzdx.begin(),drzdx.end(),0.0);
            model_drzdx(drzdx.data(),t,x.data(),udata->p(),udata->k());
        }
        
        /** model specific implementation of fdrzdx
         * param[out] drzdx partial derivative of root output rz w.r.t. model states x
         * param[in] t current time
         * param[in] x current state
         * param[in] p parameter vector
         * param[in] k constant vector
         **/
        virtual void model_drzdx(double *drzdx, const realtype t, const realtype *x, const realtype *p, const realtype *k) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** State update functions for events
         * @param[in] ie event index
         * @param[in] udata object with user input
         */
        virtual void fdeltax(const int ie, const realtype t, const AmiVector x,
                             const AmiVector xdot, const AmiVector xdot_old, const UserData *udata) {
            std::fill(deltax.begin(),deltax.end(),0.0);
            model_deltax(deltax.data(),t,x.data(),udata->p(),udata->k(),ie,xdot.data(),xdot_old.data());
        }
        
        /** model specific implementation of fdeltax
         * param[out] deltax state update
         * param[in] t current time
         * param[in] x current state
         * param[in] p parameter vector
         * param[in] k constant vector
         * param[in] ie event index
         * param[in] xdot new model right hand side
         * param[in] xdot_old previous model right hand side
         **/
        virtual void model_deltax(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k,
                                 const int ie, const realtype *xdot, const realtype *xdot_old) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** Sensitivity update functions for events, total derivative
         * @param[in] ie event index
         * @param[in] udata object with user input
         */
        virtual void fdeltasx(const int ie, const realtype t, const AmiVector x, const AmiVectorArray sx,
                              const AmiVector xdot, const AmiVector xdot_old, const UserData *udata) {
            std::fill(deltasx.begin(),deltasx.end(),0.0);
            for(int ip = 0; ip < udata->nplist(); ip++)
                model_deltasx(deltasx.data(),t,x.data(),udata->p(),udata->k(),
                             udata->plist(ip),ie,xdot.data(),xdot_old.data(),sx.data(ip),stau.data());
        }
        
        /** model specific implementation of fdeltasx
         * param[out] deltasx sensitivity update
         * param[in] t current time
         * param[in] x current state
         * param[in] p parameter vector
         * param[in] k constant vector
         * param[in] ip sensitivity index
         * param[in] ie event index
         * param[in] xdot new model right hand side
         * param[in] xdot_old previous model right hand side
         * param[in] sx state sensitivity
         * param[in] stau event-time sensitivity
         **/
        virtual void model_deltasx(double *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k,
                                   const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx,
                                   const realtype *stau) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** Adjoint state update functions for events
         * @param[in] ie event index
         * @param[in] udata object with user input
         */
        virtual void fdeltaxB(const int ie, const realtype t, const AmiVector x, const AmiVector xB,
                              const AmiVector xdot, const AmiVector xdot_old, const UserData *udata) {
            std::fill(deltaxB.begin(),deltaxB.end(),0.0);
            model_deltaxB(deltaxB.data(),t,x.data(),udata->p(),udata->k(),ie,xdot.data(),xdot_old.data(),xB.data());
        }
        
        /** model specific implementation of fdeltaxB
         * param[out] deltaxB adjoint state update
         * param[in] t current time
         * param[in] x current state
         * param[in] p parameter vector
         * param[in] k constant vector
         * param[in] ie event index
         * param[in] xdot new model right hand side
         * param[in] xdot_old previous model right hand side
         * param[in] xB current adjoint state
         **/
        virtual void model_deltaxB(double *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k,
                                  const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** Quadrature state update functions for events
         * @param[in] ie event index
         * @param[in] udata object with user input
         */
        virtual void fdeltaqB(const int ie, const realtype t, const AmiVector x, const AmiVector xB,
                              const AmiVector xdot, const AmiVector xdot_old, const AmiVector qBdot, const UserData *udata) {
            std::fill(deltaqB.begin(),deltaqB.end(),0.0);
            for(int ip = 0; ip < udata->nplist(); ip++)
                model_deltaqB(deltaqB.data(),t,x.data(),udata->p(),udata->k(),
                              udata->plist(ip),ie,xdot.data(),xdot_old.data(),xB.data(),qBdot.data());
        }
        
        /** model specific implementation of fdeltasx
         * param[out] deltasx sensitivity update
         * param[in] t current time
         * param[in] x current state
         * param[in] p parameter vector
         * param[in] k constant vector
         * param[in] ip sensitivity index
         * param[in] ie event index
         * param[in] xdot new model right hand side
         * param[in] xdot_old previous model right hand side
         * param[in] xB adjoint state
         * param[in] qBdot right hand side of quadradature
         **/
        virtual void model_deltaqB(double *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k,
                                   const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB,
                                   const realtype *qBdot) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        /** Standard deviation of measurements
         * @param[in] udata object with user input
         */
        virtual void fsigma_y(const realtype t, const UserData *udata) {
            std::fill(sigmay.begin(),sigmay.end(),0.0);
            model_sigma_y(sigmay.data(),t,udata->p(),udata->k());
        }
        
        /** model specific implementation of fsigmay
         * param[out] sigmay standard deviation of measurements
         * param[in] t current time
         * param[in] p parameter vector
         * param[in] k constant vector
         **/
        virtual void model_sigma_y(double *sigmay, const realtype t, const realtype *p, const realtype *k) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** partial derivative of standard deviation of measurements w.r.t. model
         * @param[in] udata object with user input
         */
        virtual void fdsigma_ydp(const realtype t, const UserData *udata) {
            std::fill(dsigmaydp.begin(),dsigmaydp.end(),0.0);
            for(int ip = 0; ip < udata->nplist(); ip++)
                model_dsigma_ydp(dsigmaydp.data(),t,udata->p(),udata->k(),udata->plist(ip));
        }
        
        /** model specific implementation of fsigmay
         * param[out] dsigmaydp partial derivative of standard deviation of measurements
         * param[in] t current time
         * param[in] p parameter vector
         * param[in] k constant vector
         **/
        virtual void model_dsigma_ydp(double *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** Standard deviation of events
         * @param[in] udata object with user input
         */
        virtual void fsigma_z(const realtype t, const UserData *udata) {
            std::fill(sigmaz.begin(),sigmaz.end(),0.0);
            model_sigma_z(sigmaz.data(),t,udata->p(),udata->k());
        }
        
        /** model specific implementation of fsigmaz
         * param[out] sigmaz standard deviation of event measurements
         * param[in] t current time
         * param[in] p parameter vector
         * param[in] k constant vector
         **/
        virtual void model_sigma_z(double *sigmaz, const realtype t, const realtype *p, const realtype *k) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** Sensitivity of standard deviation of events measurements w.r.t. model parameters p
         * @param[in] udata object with user input
         */
        virtual void fdsigma_zdp(const realtype t, const UserData *udata) {
            std::fill(dsigmazdp.begin(),dsigmazdp.end(),0.0);
            for(int ip = 0; ip < udata->nplist(); ip++)
                model_dsigma_zdp(dsigmazdp.data(),t,udata->p(),udata->k(),udata->plist(ip));
        }
        
        /** model specific implementation of fsigmaz
         * param[out] dsigmazdp partial derivative of standard deviation of event measurements
         * param[in] t current time
         * param[in] p parameter vector
         * param[in] k constant vector
         **/
        virtual void model_dsigma_zdp(double *dsigmazdp, const realtype t, const realtype *p, const realtype *k, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** negative log-likelihood of measurements y
         * @param[in] it timepoint index
         * @param[in] udata object with user input
         * @param[in] rdata pointer to return data object
         * @param[in] edata pointer to experimental data object
         */
        virtual void fJy(std::vector<double> Jy,const int it, const ReturnData *rdata, const UserData *udata, const ExpData *edata) {
            std::vector<double> nllh(nJ,0.0);
            std::vector<double> y = gety(it,rdata,udata);
            std::vector<double> my = getmy(it,edata,udata);
            for(int iytrue = 0; iytrue < nytrue; iytrue++){
                if(!amiIsNaN(edata->my[iytrue* udata->nt()+it])){
                    std::fill(nllh.begin(),nllh.end(),0.0);
                    model_Jy(nllh.data(),udata->p(),udata->k(),y.data(),sigmay.data(),my.data());
                    for(int iJ = 0; iJ < nJ; iJ++){
                        Jy.at(iJ) += nllh.at(iJ);
                    }
                }
            }
        }
        
        /** model specific implementation of fJy
         * param[out] nllh negative log-likelihood for measurements y
         * param[in] p parameter vector
         * param[in] k constant vector
         * param[in] y model output at timepoint
         * param[in] sigmay measurement standard deviation at timepoint
         * param[in] my measurements at timepoint
         **/
        virtual void model_Jy(double *nllh, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** negative log-likelihood of event-resolved measurements z
         * @param[in] nroots event index
         * @param[in] udata object with user input
         * @param[in] rdata pointer to return data object
         * @param[in] edata pointer to experimental data object
         */
        virtual void fJz(std::vector<double> Jz, const int nroots, const ReturnData *rdata, const UserData *udata, const ExpData *edata) {
            std::vector<double> nllh(nJ,0.0);
            std::vector<double> z = getz(nroots,rdata,udata);
            std::vector<double> mz = getmz(nroots,edata,udata);
            for(int iztrue = 0; iztrue < nztrue; iztrue++){
                if(!amiIsNaN(edata->mz[nroots + udata->nme()*iztrue])){
                    std::fill(nllh.begin(),nllh.end(),0.0);
                    model_Jz(nllh.data(),udata->p(),udata->k(),z.data(),sigmaz.data(),mz.data());
                    for(int iJ = 0; iJ < nJ; iJ++){
                        Jz.at(iJ) += nllh.at(iJ);
                    }
                }
            }
        }
        
        /** model specific implementation of fJz
         * param[out] nllh negative log-likelihood for event measurements z
         * param[in] p parameter vector
         * param[in] k constant vector
         * param[in] z model event output at timepoint
         * param[in] sigmaz event measurement standard deviation at timepoint
         * param[in] mz event measurements at timepoint
         **/
        virtual void model_Jz(double *nllh, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** regularization of negative log-likelihood with roots of event-resolved
         * measurements rz
         * @param[in] nroots event index
         * @param[in] udata object with user input
         * @param[in] rdata pointer to return data object
         * @param[in] edata pointer to experimental data object
         */
        virtual void fJrz(std::vector<double> Jz, const int nroots, const ReturnData *rdata, const UserData *udata, const ExpData *edata) {
            std::vector<double> nllh(nJ,0.0);
            std::vector<double> rz = getrz(nroots,rdata,udata);
            for(int iztrue = 0; iztrue < nztrue; iztrue++){
                if(!amiIsNaN(edata->mz[nroots + udata->nme()*iztrue])){
                    std::fill(nllh.begin(),nllh.end(),0.0);
                    model_Jrz(nllh.data(),udata->p(),udata->k(),rz.data(),sigmaz.data());
                    for(int iJ = 0; iJ < nJ; iJ++){
                        Jz.at(iJ) += nllh.at(iJ);
                    }
                }
            }
        }
        
        /** model specific implementation of fJrz
         * param[out] nllh regularization for event measurements z
         * param[in] p parameter vector
         * param[in] k constant vector
         * param[in] z model event output at timepoint
         * param[in] sigmaz event measurement standard deviation at timepoint
         **/
        virtual void model_Jrz(double *nllh, const realtype *p, const realtype *k, const double *z, const double *sigmaz) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** partial derivative of time-resolved measurement negative log-likelihood Jy
         * @param[in] it timepoint index
         * @param[in] udata object with user input
         * @param[in] rdata pointer to return data object
         * @param[in] edata pointer to experimental data object
         */
        virtual void fdJydy(const int it, const ReturnData *rdata,
                            const UserData *udata, const ExpData *edata) {
            std::vector<double> y = gety(it,rdata,udata);
            std::vector<double> my = getmy(it,edata,udata);
            std::vector<double> dJydy_slice(nJ*nytrue, 0.0);
            std::fill(dJydy.begin(),dJydy.end(),0.0);
            for(int iytrue = 0; iytrue < nytrue; iytrue++){
                if(!amiIsNaN(edata->my[iytrue* udata->nt()+it])){
                    std::fill(dJydy_slice.begin(),dJydy_slice.end(),0.0);
                    model_dJydy(dJydy_slice.data(),udata->p(),udata->k(),y.data(),sigmay.data(),my.data());
                    // TODO: fix slicing here such that slicing is no longer necessary in sJy
                    for(int iJ = 0; iJ < nJ; iJ++){
                        for(int iy = 0; iy < ny; iy++ )
                        dJydy.at(iytrue+(iJ+iy*nJ)) = dJydy_slice.at(iJ+iy*nJ);
                    }
                }
            }
        }
        
        /** model specific implementation of fdJydy
         * param[out] dJydy partial derivative of time-resolved measurement negative log-likelihood Jy
         * param[in] p parameter vector
         * param[in] k constant vector
         * param[in] y model output at timepoint
         * param[in] sigmay measurement standard deviation at timepoint
         * param[in] my measurement at timepoint
         **/
        virtual void model_dJydy(double *dJydy, const realtype *p, const realtype *k,
                                 const double *y, const double *sigmay, const double *my) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** Sensitivity of time-resolved measurement negative log-likelihood Jy
         * w.r.t. standard deviation sigma
         * @param[in] it timepoint index
         * @param[in] udata object with user input
         * @param[in] rdata pointer to return data object
         * @param[in] edata pointer to experimental data object
         */
        virtual void fdJydsigma(const int it, const ReturnData *rdata,
                                const UserData *udata, const ExpData *edata) {
            std::vector<double> y = gety(it,rdata,udata);
            std::vector<double> my = getmy(it,edata,udata);
            std::vector<double> dJydsigma_slice(nJ*nytrue, 0.0);
            std::fill(dJydsigma.begin(),dJydsigma.end(),0.0);
            for(int iytrue = 0; iytrue < nytrue; iytrue++){
                if(!amiIsNaN(edata->my[iytrue* udata->nt()+it])){
                    std::fill(dJydsigma_slice.begin(),dJydsigma_slice.end(),0.0);
                    model_dJydsigma(dJydsigma_slice.data(),udata->p(),udata->k(),y.data(),sigmay.data(),my.data());
                    // TODO: fix slicing here such that slicing is no longer necessary in sJy
                    for(int iJ = 0; iJ < nJ; iJ++){
                        for(int iy = 0; iy < ny; iy++ )
                            dJydsigma.at(iytrue+(iJ+iy*nJ)) = dJydsigma_slice.at(iJ+iy*nJ);
                    }
                }
            }
        }
        
        /** model specific implementation of fdJydsigma
         * param[out] dJydsigma Sensitivity of time-resolved measurement 
         * negative log-likelihood Jy w.r.t. standard deviation sigmay
         * param[in] p parameter vector
         * param[in] k constant vector
         * param[in] y model output at timepoint
         * param[in] sigmaz measurement standard deviation at timepoint
         * param[in] my measurement at timepoint
         **/
        virtual void model_dJydsigma(double *dJydsigma, const realtype *p, const realtype *k,
                                 const double *y, const double *sigmay, const double *my) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** partial derivative of event measurement negative log-likelihood Jz
         * @param[in] nroots event index
         * @param[in] udata object with user input
         * @param[in] rdata pointer to return data object
         * @param[in] edata pointer to experimental data object
         */
        virtual void fdJzdz(const int nroots, const ReturnData *rdata,
                            const UserData *udata, const ExpData *edata) {
            std::vector<double> z = getz(nroots,rdata,udata);
            std::vector<double> mz = getmz(nroots,edata,udata);
            std::vector<double> dJzdz_slice(nJ*nztrue, 0.0);
            std::fill(dJzdz.begin(),dJzdz.end(),0.0);
            for(int iztrue = 0; iztrue < nztrue; iztrue++){
                if(!amiIsNaN(edata->mz[iztrue*udata->nme()+nroots])){
                    std::fill(dJzdz_slice.begin(),dJzdz_slice.end(),0.0);
                    model_dJzdz(dJzdz_slice.data(),udata->p(),udata->k(),z.data(),sigmaz.data(),mz.data());
                    // TODO: fix slicing here such that slicing is no longer necessary in sJz
                    for(int iJ = 0; iJ < nJ; iJ++){
                        for(int iz = 0; iz < nz; iz++ )
                            dJzdz.at(iztrue+(iJ+iz*nJ)) = dJzdz_slice.at(iJ+iz*nJ);
                    }
                }
            }
        }
        
        /** model specific implementation of fdJzdz
         * param[out] dJzdz partial derivative of event measurement negative log-likelihood Jz
         * param[in] p parameter vector
         * param[in] k constant vector
         * param[in] z model event output at timepoint
         * param[in] sigmaz event measurement standard deviation at timepoint
         * param[in] mz event measurement at timepoint
         **/
        virtual void model_dJzdz(double *dJzdz, const realtype *p, const realtype *k,
                                 const double *z, const double *sigmaz, const double *mz) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** Sensitivity of event measurement negative log-likelihood Jz
         * w.r.t. standard deviation sigmaz
         * @param[in] nroots event index
         * @param[in] udata object with user input
         * @param[in] rdata pointer to return data object
         * @param[in] edata pointer to experimental data object
         */
        virtual void fdJzdsigma(const int nroots, const ReturnData *rdata,
                                const UserData *udata, const ExpData *edata) {
            std::vector<double> z = getz(nroots,rdata,udata);
            std::vector<double> mz = getmz(nroots,edata,udata);
            std::vector<double> dJzdsigma_slice(nJ*nztrue, 0.0);
            std::fill(dJzdsigma.begin(),dJzdsigma.end(),0.0);
            for(int iztrue = 0; iztrue < nztrue; iztrue++){
                if(!amiIsNaN(edata->mz[iztrue*udata->nme()+nroots])){
                    std::fill(dJzdsigma_slice.begin(),dJzdsigma_slice.end(),0.0);
                    model_dJzdsigma(dJzdsigma_slice.data(),udata->p(),udata->k(),z.data(),sigmaz.data(),mz.data());
                    // TODO: fix slicing here such that slicing is no longer necessary in sJz
                    for(int iJ = 0; iJ < nJ; iJ++){
                        for(int iz = 0; iz < nz; iz++ )
                            dJzdsigma.at(iztrue+(iJ+iz*nJ)) = dJzdsigma_slice.at(iJ+iz*nJ);
                    }
                }
            }
        }
        
        /** model specific implementation of fdJzdsigma
         * param[out] dJzdsigma Sensitivity of event measurement
         * negative log-likelihood Jz w.r.t. standard deviation sigmaz
         * param[in] p parameter vector
         * param[in] k constant vector
         * param[in] z model event output at timepoint
         * param[in] sigmaz event measurement standard deviation at timepoint
         * param[in] mz event measurement at timepoint
         **/
        virtual void model_dJzdsigma(double *dJzdsigma, const realtype *p, const realtype *k,
                                     const double *z, const double *sigmaz, const double *mz) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** partial derivative of event measurement negative log-likelihood Jz
         * @param[in] nroots event index
         * @param[in] udata object with user input
         * @param[in] rdata pointer to return data object
         * @param[in] edata pointer to experimental data object
         */
        virtual void fdJrzdz(const int nroots, const ReturnData *rdata,
                            const UserData *udata, const ExpData *edata) {
            std::vector<double> rz = getz(nroots,rdata,udata);
            std::vector<double> mz = getmz(nroots,edata,udata);
            std::vector<double> dJrzdz_slice(nJ*nztrue, 0.0);
            std::fill(dJrzdz.begin(),dJrzdz.end(),0.0);
            for(int iztrue = 0; iztrue < nztrue; iztrue++){
                if(!amiIsNaN(edata->mz[iztrue*udata->nme()+nroots])){
                    std::fill(dJrzdz_slice.begin(),dJrzdz_slice.end(),0.0);
                    model_dJzdz(dJrzdz_slice.data(),udata->p(),udata->k(),rz.data(),sigmaz.data());
                    // TODO: fix slicing here such that slicing is no longer necessary in sJz
                    for(int iJ = 0; iJ < nJ; iJ++){
                        for(int iz = 0; iz < nz; iz++)
                            dJrzdz.at(iztrue+(iJ+iz*nJ)) = dJrzdz_slice.at(iJ+iz*nJ);
                    }
                }
            }
        }
        
        /** model specific implementation of fdJrzdz
         * param[out] dJrzdz partial derivative of event penalization Jrz
         * param[in] p parameter vector
         * param[in] k constant vector
         * param[in] rz model root output at timepoint
         * param[in] sigmaz event measurement standard deviation at timepoint
         **/
        virtual void model_dJrzdz(double *dJrzdz, const realtype *p, const realtype *k,
                                 const double *rz, const double *sigmaz) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** Sensitivity of event measurement negative log-likelihood Jz
         * w.r.t. standard deviation sigmaz
         * @param[in] nroots event index
         * @param[in] udata object with user input
         * @param[in] rdata pointer to return data object
         * @param[in] edata pointer to experimental data object
         */
        virtual void fdJrzdsigma(const int nroots,const ReturnData *rdata,
                                const UserData *udata, const ExpData *edata) {
            std::vector<double> rz = getrz(nroots,rdata,udata);
            std::vector<double> dJrzdsigma_slice(nJ*nztrue, 0.0);
            std::fill(dJrzdsigma.begin(),dJrzdsigma.end(),0.0);
            for(int iztrue = 0; iztrue < nztrue; iztrue++){
                if(!amiIsNaN(edata->mz[iztrue*udata->nme()+nroots])){
                    std::fill(dJrzdsigma_slice.begin(),dJrzdsigma_slice.end(),0.0);
                    model_dJrzdsigma(dJrzdsigma_slice.data(),udata->p(),udata->k(),rz.data(),sigmaz.data());
                    // TODO: fix slicing here such that slicing is no longer necessary in sJz
                    for(int iJ = 0; iJ < nJ; iJ++){
                        for(int iz = 0; iz < nz; iz++ )
                            dJrzdsigma.at(iztrue+(iJ+iz*nJ)) = dJrzdsigma_slice.at(iJ+iz*nJ);
                    }
                }
            }
        }
        
        /** model specific implementation of fdJrzdsigma
         * param[out] dJzdsigma Sensitivity of event penalization Jz w.r.t.
         * standard deviation sigmaz
         * param[in] p parameter vector
         * param[in] k constant vector
         * param[in] rz model root output at timepoint
         * param[in] sigmaz event measurement standard deviation at timepoint
         **/
        virtual void model_dJrzdsigma(double *dJrzdsigma, const realtype *p, const realtype *k,
                                     const double *rz, const double *sigmaz) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        virtual ~Model() {
            SparseDestroyMat(J);
        };
        
        // Generic implementations
        void fsy(const int it, const TempData *tdata, ReturnData *rdata);
        
        void fsz_tf(const int ie, const TempData *tdata, ReturnData *rdata);
        
        void fsJy(const int it, const TempData *tdata, ReturnData *rdata);
        
        void fdJydp(const int it, TempData *tdata, const ExpData *edata,
                   const ReturnData *rdata);
        
        void fdJydx(const int it, TempData *tdata, const ExpData *edata);
        
        void fsJz(const int ie, TempData *tdata, const ReturnData *rdata);
        
        void fdJzdp(const int ie, TempData *tdata, const ExpData *edata,
                   const ReturnData *rdata);
        
        void fdJzdx(const int ie, TempData *tdata, const ExpData *edata);
        
        void initialize(const UserData *udata, TempData *tdata);
        
        void initializeStates(const double *x0data, TempData *tdata);
        
        void initHeaviside(TempData *tdata);
        
    protected:
        int checkVals(const int N,const realtype *array, const char* fun){
            for(int idx = 0; idx < N; idx++) {
                if(amiIsNaN(array[idx])) {
                    warnMsgIdAndTxt("AMICI:mex:fJDiag:NaN","AMICI replaced a NaN value at index (%i) of (%i) in (%s)! Aborting simulation ... ",idx,N,fun);
                    return(AMICI_ERROR);
                }
                if(amiIsInf(array[idx])) {
                    warnMsgIdAndTxt("AMICI:mex:fJDiag:Inf","AMICI encountered an Inf value at index (%i) of (%i) in (%s)! Aborting simulation ... ",idx,N,fun);
                    return(AMICI_ERROR);
                }
            }
            return(AMICI_SUCCESS);
        }
        
        /**
         * @brief Recurring terms in xdot
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] udata object with user input
         */
        virtual void fw(const realtype t, const N_Vector x, const UserData *udata) {
            std::fill(w.begin(),w.end(),0.0);
            model_w(w.data(),t,N_VGetArrayPointer(x),udata->p(),udata->k());
        }
        
        /** model specific implementation of fw
         * @param[out] w Recurring terms in xdot
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] p parameter vector
         * @param[in] k constants vector
         */
        virtual void model_w(realtype *w, const realtype t, const realtype *x, const realtype *p,
                                const realtype *k) {};
        
        /**
         * @brief Recurring terms in xdot, parameter derivative
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] udata object with user input
         */
        virtual void fdwdp(const realtype t, const N_Vector x, const UserData *udata) {
            fw(t,x,udata);
            std::fill(dwdp.begin(),dwdp.end(),0.0);
            model_dwdp(dwdp.data(),t,N_VGetArrayPointer(x),udata->p(),udata->k(),w.data());
        }
        
        /** model specific implementation of dwdp
         * @param[out] dwdp Recurring terms in xdot, parameter derivative
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] p parameter vector
         * @param[in] k constants vector
         * @param[in] w vector with helper variables
         */
        virtual void model_dwdp(realtype *dwdp, const realtype t, const realtype *x, const realtype *p,
                                const realtype *k, const realtype *w) {};
        
        /**
         * @brief Recurring terms in xdot, state derivative
         * @param[in] t timepoint
         * @param[in] x Vector with the states @type N_Vector
         * @param[in] udata object with user input
         */
        virtual void fdwdx(const realtype t, const N_Vector x, const UserData *udata) {
            fw(t,x,udata);
            std::fill(dwdx.begin(),dwdx.end(),0.0);
            model_dwdx(dwdx.data(),t,N_VGetArrayPointer(x),udata->p(),udata->k(),w.data());
        }
        
        /** model specific implementation of dwdx
         * @param[out] dwdx Recurring terms in xdot, state derivative
         * @param[in] t timepoint
         * @param[in] x Vector with the states
         * @param[in] p parameter vector
         * @param[in] k constants vector
         * @param[in] w vector with helper variables
         */
        virtual void model_dwdx(realtype *dwdx, const realtype t, const realtype *x, const realtype *p,
                                const realtype *k, const realtype *w) {};
        
        /* Model dimensions */
        int nplist;
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
        const std::vector<int> idlist;
        
    private:
        
        /** create my slice at timepoint
         * @param[in] it timepoint index
         * @param[in] udata object with user input
         * @param[in] edata pointer to experimental data object
         */
        std::vector<double> getmy(const int it, const ExpData *edata, const UserData *udata) {
            std::vector<double> my(nytrue, 0.0);
            for(int iytrue = 0; iytrue < nytrue; iytrue++){
                my.at(iytrue) = edata->my[it + udata->nt()*iytrue];
            }
            return my;
        }
        
        /** create y slice at timepoint
         * @param[in] it timepoint index
         * @param[in] udata object with user input
         * @param[in] rdata pointer to return data object
         */
        std::vector<double> gety(const int it, const ReturnData *rdata, const UserData *udata) {
            std::vector<double> y(ny, 0.0);
            for(int iy = 0; iy < ny; iy++){
                y.at(iy) = rdata->y[it + udata->nt()*iy];
            }
            return y;
        }
        
        /** create mz slice at event
         * @param[in] nroots event occurence
         * @param[in] udata object with user input
         * @param[in] edata pointer to experimental data object
         */
        std::vector<double> getmz(const int nroots, const ExpData *edata, const UserData *udata) {
            std::vector<double> mz(nztrue, 0.0);
            for(int iztrue = 0; iztrue < nztrue; iztrue++){
                mz.at(iztrue) = edata->mz[nroots + udata->nme()*iztrue];
            }
            return mz;
        }
        
        /** create z slice at event
         * @param[in] nroots event occurence
         * @param[in] udata object with user input
         * @param[in] rdata pointer to return data object
         */
        std::vector<double> getz(const int nroots, const ReturnData *rdata, const UserData *udata) {
            std::vector<double> z(nz, 0.0);
            for(int iz = 0; iz < nz; iz++){
                z.at(iz) = rdata->z[nroots+udata->nme()*iz];
            }
            return z;
        }
        
        /** create rz slice at event
         * @param[in] nroots event occurence
         * @param[in] udata object with user input
         * @param[in] rdata pointer to return data object
         */
        std::vector<double> getrz(const int nroots, const ReturnData *rdata, const UserData *udata) {
            std::vector<double> rz(nz, 0.0);
            for(int iz = 0; iz < nz; iz++){
                rz.at(iz) = rdata->rz[nroots+udata->nme()*iz]
            }
            return rz;
        }

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
        SlsMat J;
        
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
        /** data standard deviation */
        std::vector<double> sigmay;
        /** parameter derivative of data standard deviation */
        std::vector<double> dsigmaydp;
        /** event standard deviation */
        std::vector<double> sigmaz;
        /** parameter derivative of event standard deviation */
        std::vector<double> dsigmazdp;
        /** tempory storage of dxdotdp data across functions */
        std::vector<realtype> dxdotdp;;
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
        /** change in x */
        std::vector<realtype> deltax;
        /** change in sx */
        std::vector<realtype> deltasx;
        /** change in xB */
        std::vector<realtype> deltaxB;
        /** change in qB */
        std::vector<realtype> deltaqB;
    };
    
} // namespace amici

#endif // MODEL_H
