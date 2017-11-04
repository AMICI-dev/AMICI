#ifndef AMICI_MODEL_H
#define AMICI_MODEL_H

#include <include/amici.h>
#include <include/amici_exception.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_direct.h>
#include <sundials/sundials_sparse.h>

#include <vector>

namespace amici {
    
    class UserData;
    class ExpData;
    
    /**
     * @brief The Model class represents an AMICI ODE model.
     * The model does not contain any data, but represents the state
     * of the model at a specific time t. The states must not always be
     * in sync, but may be updated asynchroneously. 
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
        nnz(nnz), nJ(nJ), ubw(ubw), lbw(lbw), o2mode(o2mode) {
        
            dJydyTmp.resize(nJ * ny, 0);
            dJydxTmp.resize(nJ * nx, 0);
            dJydsigmaTmp.resize(nJ * ny, 0);
            dJzdzTmp.resize(nJ * nz, 0);
            dJzdxTmp.resize(nJ * nx, 0);
            dJzdsigmaTmp.resize(nJ * nz, 0);
            dJrzdsigmaTmp.resize(nJ * nz, 0);
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
        virtual void fx0(UserData *udata) {
            x0.reset()
            model_x0(x0.data,udata->tstart,udata.p(),udata.k());
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
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
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
        virtual void fsx0(const UserData *udata) {
            sx0.reset();
            for(ip = 0; ip<udata->nplist; ip++)
                model_sx0(sx0.data(ip),udata->tstart,x.data(),udata.p(),udata.k(),udata->plist[ip]);
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
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
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
        virtual void fstau(const int ie, const UserData *udata) {
            std::fill(stau.begin(),stau.end(),0.0);
            for(int ip = 0; ip < udata->nplist; ip++){
                model_stau(stau.data(),t,x.data(),udata->p(),udata->k(),sx.data(ip),udata->plist[ip],ie)
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
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        
        /** Observables / measurements
         * @param[in] it timepoint index
         * @param[out] rdata pointer to return data object
         * @param[in] udata object with user input
         */
        virtual void fy(int it, ReturnData *rdata
                        const UserData *udata) {
            std::vector<double> yreturn(ny,0.0);
            model_y(yreturn,t,x,p,k);
            for(int iy; iy < ny; iy++)
                rdata->y[it + udata->nt*iy] = yreturn.at(iy);
        }
        
        /** model specific implementation of fy
         * param[out] y model output at current timepoint
         * param[in] t current time
         * param[in] x current state
         * param[in] p parameter vector
         * param[in] k constant vector
         **/
        virtual void model_y(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k) {
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** partial derivative of observables y w.r.t. model parameters p
         * @param[in] udata object with user input
         */
        virtual void fdydp(const UserData *udata) {
            std::fill(dydp.begin(),dydp.end(),0.0);
            for(int ip = 0; ip < udata->nplist; ip++){
                model_dydp(dydp.data(),t,x.data(),udata->p,udata->k,udata->plist[ip])
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
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** partial derivative of observables y w.r.t. state variables x
         const UserData *udata
         */
        virtual void fdydx(const UserData *udata) {
            std::fill(dydx.begin(),dydx.end(),0.0);
            model_dypdx(dydx.data(),t,x.data(),udata->p,udata->k)
        }
        
        /** model specific implementation of fdydx
         * param[out] dydx partial derivative of observables y w.r.t. model states x
         * param[in] t current time
         * param[in] x current state
         * param[in] p parameter vector
         * param[in] k constant vector
         **/
        virtual void model_dydx(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k) {
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** Event-resolved output
         * @param[in] nroots number of events for event index
         * @param[out] rdata pointer to return data object
         * @param[in] udata object with user input
         */
        virtual void fz(const int nroots, ReturnData *rdata,
                        const UserData *udata) {
            std::vector<double> zreturn(nz,0.0);
            model_z(zreturn.data(),t,x.data(),udata->p(),udata->k());
            for(int iz; iz < nz; iz++) {
                rdata->z[nroots+udata->nmaxevent*iz] = rzreturn.at(iz);
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
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** Sensitivity of z, total derivative
         * @param[in] nroots number of events for event index
         * @param[out] rdata pointer to return data object
         * @param[in] udata object with user input
         */
        virtual void fsz(const int nroots, ReturnData *rdata
                         const UserData *udata) {
            for(int ip; ip < udata->nplist;  ip++ ){
                std::vector<double> szreturn(nz,0.0);
                model_sz(szreturn.data(),t,x.data(),udata->p(),udata->k(),sx.data(ip),udata->plist[ip]);
                for(int iz; iz < nz; iz++) {
                    rdata->sz[tdata->nroots[ie]+udata->nmaxevent*(ip*model->nz + iz)] = szreturn.at(iz);
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
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** Event root function of events (equal to froot but does not include
         * non-output events)
         * @param[in] nroots number of events for event index
         * @param[in] udata object with user input
         * @param[out] rdata pointer to return data object
         */
        virtual void frz(const int nroots, ReturnData *rdata
                         const UserData *udata) {
            std::vector<double> rzreturn(nz,0.0);
            model_rz(rzreturn.data(),t,x.data(),udata->p(),udata->k());
            for(int iz; iz < nz; iz++) {
                rdata->rz[nroots+udata->nmaxevent*iz] = rzreturn.at(iz);
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
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** Sensitivity of rz, total derivative
         * @param[in] nroots number of events for event index
         * @param[out] rdata pointer to return data object
         * @param[in] udata object with user input
         */
        virtual void fsrz(const int nroots, ReturnData *rdata
                          const UserData *udata) {
            for(int ip; ip < udata->nplist;  ip++ ){
                std::vector<double> srzreturn(nz,0.0);
                model_srz(srzreturn.data(),t,x.data(),udata->p(),udata->k(),sx.data(ip),udata->plist[ip]);
                for(int iz; iz < nz; iz++) {
                    rdata->srz[tdata->nroots[ie]+udata->nmaxevent*(ip*model->nz + iz)] = srzreturn.at(iz);
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
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** partial derivative of event-resolved output z w.r.t. to model parameters p
         * @param[in] udata object with user input
         */
        virtual void fdzdp(const UserData *udata) {
            std::fill(dxdp.begin(),dzdp.end(),0.0);
            for(int ip = 0; ip < udata->nplist; ip++){
                model_dzdp(dzdp.data(),t,x.data(),udata->p,udata->k,udata->plist[ip]);
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
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** partial derivative of event-resolved output z w.r.t. to model states x
         * @param[in] udata object with user input
         */
        virtual void fdzdx(const UserData *udata) {
            std::fill(dzdx.begin(),dzdx.end(),0.0);
            model_dzdx(dzdx.data(),t,x.data(),udata->p,udata->k);
        }
        
        /** model specific implementation of fdzdx
         * param[out] dzdx partial derivative of event-resolved output z w.r.t. model states x
         * param[in] t current time
         * param[in] x current state
         * param[in] p parameter vector
         * param[in] k constant vector
         **/
        virtual void model_dzdx(double *dzdx, const realtype t, const realtype *x, const realtype *p, const realtype *k) {
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** Sensitivity of event-resolved root output w.r.t. to model parameters p
         * @param[in] udata object with user input
         */
        virtual void fdrzdp(const UserData *udata) {
            std::fill(drzdp.begin(),drzdp.end(),0.0);
            for(int ip = 0; ip < udata->nplist; ip++){
                model_drzdp(drzdp.data(),t,x.data(),udata->p,udata->k,udata->plist[ip]);
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
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** Sensitivity of event-resolved measurements rz w.r.t. to model states x
         * @param[in] udata object with user input
         */
        virtual void fdrzdx(const UserData *udata) {
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
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** State update functions for events
         * @param[in] ie event index
         * @param[in] udata object with user input
         */
        virtual void fdeltax(const int ie, const UserData *udata) {
            std::fill(deltax.begin(),deltax.end(),0.0)
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
        virtual void model_deltax(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k
                                 const int ie, const realtype *xdot, const realtype *xdot_old) {
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** Sensitivity update functions for events, total derivative
         * @param[in] ie event index
         * @param[in] udata object with user input
         */
        virtual void fdeltasx(const int ie, const UserData *udata) {
            std::fill(deltasx.begin(),deltasx.end(),0.0)
            for(int ip = 0; ip < udata->nplist; ip++)
                model_deltasx(deltasx.data(),t,x.data(),udata->p(),udata->k(),
                             udata->plist[ip],ie,xdot.data(),xdot_old.data(),sx.data(ip),stau.data());
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
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** Adjoint state update functions for events
         * @param[in] ie event index
         * @param[in] udata object with user input
         */
        virtual void fdeltaxB(const int ie, const UserData *udata) {
            std::fill(deltaxB.begin(),deltaxB.end(),0.0)
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
        virtual void model_deltaxB(double *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k
                                  const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB) {
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** Quadrature state update functions for events
         * @param[in] ie event index
         * @param[in] udata object with user input
         */
        virtual void fdeltaqB(const int ie, const UserData *udata) {
            std::fill(deltaqB.begin(),deltaqB.end(),0.0)
            for(int ip = 0; ip < udata->nplist; ip++)
                model_deltaqB(deltaqB.data(),t,x.data(),udata->p(),udata->k(),
                              udata->plist[ip],ie,xdot.data(),xdot_old.data(),xB.data(),qBdot.data());
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
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        /** Standard deviation of measurements
         * @param[in] udata object with user input
         */
        virtual void fsigma_y(const UserData *udata) {
            std::fill(sigmay.begin(),sigmay.end(),0.0)
            model_sigma_y(sigmay.data(),t,udata->p(),udata->k());
        }
        
        /** model specific implementation of fsigmay
         * param[out] sigmay standard deviation of measurements
         * param[in] t current time
         * param[in] p parameter vector
         * param[in] k constant vector
         **/
        virtual void model_sigma_y(double *sigmay, const realtype t, const realtype *p, const realtype *k) {
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** partial derivative of standard deviation of measurements w.r.t. model
         * @param[in] udata object with user input
         */
        virtual void fdsigma_ydp(const UserData *udata) {
            std::fill(dsigmaydp.begin(),dsigmaydp.end(),0.0)
            for(int ip = 0; ip < udata->nplist; ip++)
                model_dsigma_ydp(dsigmaydp.data(),t,udata->p(),udata->k(),udata->plist[ip]);
        }
        
        /** model specific implementation of fsigmay
         * param[out] dsigmaydp partial derivative of standard deviation of measurements
         * param[in] t current time
         * param[in] p parameter vector
         * param[in] k constant vector
         **/
        virtual void model_dsigma_ydp(double *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const int ip) {
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** Standard deviation of events
         * @param[in] udata object with user input
         */
        virtual void fsigma_z(const UserData *udata) {
            std::fill(sigmaz.begin(),sigmaz.end(),0.0)
            model_sigma_z(sigmaz.data(),t,udata->p(),udata->k());
        }
        
        /** model specific implementation of fsigmaz
         * param[out] sigmaz standard deviation of event measurements
         * param[in] t current time
         * param[in] p parameter vector
         * param[in] k constant vector
         **/
        virtual void model_sigma_z(double *sigmaz, const realtype t, const realtype *p, const realtype *k) {
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** Sensitivity of standard deviation of events measurements w.r.t. model parameters p
         * @param[in] udata object with user input
         */
        virtual void fdsigma_zdp(const UserData *udata) {
            std::fill(dsigmazdp.begin(),dsigmazdp.end(),0.0)
            for(int ip = 0; ip < udata->nplist; ip++)
                model_dsigma_zdp(dsigmazdp.data(),t,udata->p(),udata->k(),udata->plist[ip]);
        }
        
        /** model specific implementation of fsigmaz
         * param[out] dsigmazdp partial derivative of standard deviation of event measurements
         * param[in] t current time
         * param[in] p parameter vector
         * param[in] k constant vector
         **/
        virtual void model_dsigma_zdp(double *dsigmazdp, const realtype t, const realtype *p, const realtype *k, const int ip) {
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** negative log-likelihood of measurements y
         * @param[in] it timepoint index
         * @param[in] udata object with user input
         * @param[in] rdata pointer to return data object
         * @param[in] edata pointer to experimental data object
         */
        virtual void fJy(const int it,const ReturnData *rdata, const UserData *udata, const ExpData *edata) {
            std::vector<double> nllh(nJ,0.0);
            std::vector<double> y = gety(it,rdata);
            std::vector<double> my = getmy(it,edata);
            for(int iytrue = 0; iytrue < nytrue; iytrue++){
                if(!amiIsNaN(edata->my[iytrue* udata->nt+it])){
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
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** negative log-likelihood of event-resolved measurements z
         * @param[in] nroots event index
         * @param[in] udata object with user input
         * @param[in] rdata pointer to return data object
         * @param[in] edata pointer to experimental data object
         */
        virtual void fJz(const int nroots,const ReturnData *rdata, const UserData *udata, const ExpData *edata) {
            std::vector<double> nllh(nJ,0.0);
            std::vector<double> z = getz(it,rdata);
            std::vector<double> mz = getmz(it,edata);
            for(int iztrue = 0; iztrue < nytrue; iytrue++){
                if(!amiIsNaN(edata->mz[nroots + udata->nmaxevent()*iztrue])){
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
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** regularization of negative log-likelihood with roots of event-resolved
         * measurements rz
         * @param[in] nroots event index
         * @param[in] udata object with user input
         * @param[in] rdata pointer to return data object
         * @param[in] edata pointer to experimental data object
         */
        virtual void fJrz(const int nroots, const ReturnData *rdata, const UserData *udata, const ExpData *edata) {
            std::vector<double> nllh(nJ,0.0);
            std::vector<double> rz = getrz(it,rdata);
            for(int iztrue = 0; iztrue < nytrue; iytrue++){
                if(!amiIsNaN(edata->mz[nroots + udata->nmaxevent()*iztrue])){
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
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** partial derivative of time-resolved measurement negative log-likelihood Jy
         * @param[in] it timepoint index
         * @param[in] udata object with user input
         * @param[in] rdata pointer to return data object
         * @param[in] edata pointer to experimental data object
         */
        virtual void fdJydy(const int it,const ReturnData *rdata,
                            const UserData *udata, const ExpData *edata) {
            std::vector<double> y = gety(it,rdata);
            std::vector<double> my = getmy(it,edata);
            std::vector<double> dJydy_slice(nJ*nytrue, 0.0);
            std::fill(dJydy.begin(),dJydy.end(),0.0);
            for(int iytrue = 0; iytrue < nytrue; iytrue++){
                if(!amiIsNaN(edata->my[iytrue* udata->nt+it])){
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
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** Sensitivity of time-resolved measurement negative log-likelihood Jy
         * w.r.t. standard deviation sigma
         * @param[in] it timepoint index
         * @param[in] udata object with user input
         * @param[in] rdata pointer to return data object
         * @param[in] edata pointer to experimental data object
         */
        virtual void fdJydsigma(const int it,const ReturnData *rdata,
                                const UserData *udata, const ExpData *edata) {
            std::vector<double> y = gety(it,rdata);
            std::vector<double> my = getmy(it,edata);
            std::vector<double> dJydsigma_slice(nJ*nytrue, 0.0);
            std::fill(dJydsigma.begin(),dJydsigma.end(),0.0);
            for(int iytrue = 0; iytrue < nytrue; iytrue++){
                if(!amiIsNaN(edata->my[iytrue* udata->nt+it])){
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
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** partial derivative of event measurement negative log-likelihood Jz
         * @param[in] nroots event index
         * @param[in] udata object with user input
         * @param[in] rdata pointer to return data object
         * @param[in] edata pointer to experimental data object
         */
        virtual void fdJzdz(const int nroots,const ReturnData *rdata,
                            const UserData *udata, const ExpData *edata) {
            std::vector<double> z = getz(nroots,rdata);
            std::vector<double> mz = getmz(nroots,edata);
            std::vector<double> dJzdz_slice(nJ*nztrue, 0.0);
            std::fill(dJzdz.begin(),dJzdz.end(),0.0);
            for(int iztrue = 0; iztrue < nztrue; iztrue++){
                if(!amiIsNaN(edata->mz[iztrue*udata->nmaxevent()+nroots])){
                    std::fill(dJzdz_slice.begin(),dJzdz_slice.end(),0.0);
                    model_dJzdz(dJzdz_slice.data(),udata->p(),udata->k(),z.data(),sigmaz.data(),mz.data());
                    // TODO: fix slicing here such that slicing is no longer necessary in sJz
                    for(int iJ = 0; iJ < nJ; iJ++){
                        for(int iz = 0; iz < nz; iy++ )
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
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** Sensitivity of event measurement negative log-likelihood Jz
         * w.r.t. standard deviation sigmaz
         * @param[in] nroots event index
         * @param[in] udata object with user input
         * @param[in] rdata pointer to return data object
         * @param[in] edata pointer to experimental data object
         */
        virtual void fdJzdsigma(const int nroots,const ReturnData *rdata,
                                const UserData *udata, const ExpData *edata) {
            std::vector<double> z = getz(nroots,rdata);
            std::vector<double> mz = getmz(nroots,edata);
            std::vector<double> dJzdsigma_slice(nJ*nztrue, 0.0);
            std::fill(dJzdsigma.begin(),dJzdsigma.end(),0.0);
            for(int iztrue = 0; iztrue < nztrue; iztrue++){
                if(!amiIsNaN(edata->mz[iztrue*udata->nmaxevent()+nroots])){
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
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        /** partial derivative of event measurement negative log-likelihood Jz
         * @param[in] nroots event index
         * @param[in] udata object with user input
         * @param[in] rdata pointer to return data object
         * @param[in] edata pointer to experimental data object
         */
        virtual void fdJrzdz(const int nroots,const ReturnData *rdata,
                            const UserData *udata, const ExpData *edata) {
            std::vector<double> rz = getz(nroots,rdata);
            std::vector<double> mz = getmz(nroots,edata);
            std::vector<double> dJrzdz_slice(nJ*nztrue, 0.0);
            std::fill(dJrzdz.begin(),dJrzdz.end(),0.0);
            for(int iztrue = 0; iztrue < nztrue; iztrue++){
                if(!amiIsNaN(edata->mz[iztrue*udata->nmaxevent()+nroots])){
                    std::fill(dJrzdz_slice.begin(),dJrzdz_slice.end(),0.0);
                    model_dJzdz(dJrzdz_slice.data(),udata->p(),udata->k(),z.data(),sigmaz.data());
                    // TODO: fix slicing here such that slicing is no longer necessary in sJz
                    for(int iJ = 0; iJ < nJ; iJ++){
                        for(int iz = 0; iz < nz; iy++ )
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
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
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
            std::vector<double> rz = getrz(nroots,rdata);
            std::vector<double> dJzdsigma_slice(nJ*nztrue, 0.0);
            std::fill(dJrzdsigma.begin(),dJrzdsigma.end(),0.0);
            for(int iztrue = 0; iztrue < nztrue; iztrue++){
                if(!amiIsNaN(edata->mz[iztrue*udata->nmaxevent()+nroots])){
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
            throw AmiException("Requested functionality is not supported as " __func__ "is not implemented for this model!");
        }
        
        virtual ~Model();
        
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
                mz.at(iztrue) = edata->mz[nroots + udata->nmaxevent()*iztrue];
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
                z.at(iz) = rdata->z[nroots+udata->nmaxevent()*iz]
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
                rz.at(iz) = rdata->rz[nroots+udata->nmaxevent()*iz]
            }
            return rz;
        }
        
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
        const int *z2event = nullptr;
        /** flag array for DAE equations */
        const int *idlist = nullptr;

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
        
        /** state vector */
        AmiVector x;
        /** old state vector */
        AmiVector x_old;
        /** array of state vectors at discontinuities*/
        AmiVector *x_disc;
        /** array of differential state vectors at discontinuities*/
        AmiVector *xdot_disc;
        /** array of old differential state vectors at discontinuities*/
        AmiVector *xdot_old_disc;
        /** differential state vector */
        AmiVector dx;
        /** old differential state vector */
        AmiVector dx_old;
        /** time derivative state vector */
        AmiVector xdot;
        /** old time derivative state vector */
        AmiVector xdot_old;
        /** adjoint state vector */
        AmiVector xB;
        /** old adjoint state vector */
        AmiVector xB_old;
        /** differential adjoint state vector */
        AmiVector dxB;
        /** quadrature state vector */
        AmiVector xQB;
        /** old quadrature state vector */
        AmiVector xQB_old;
        /** sensitivity state vector array */
        AmiVectorArray *sx;
        /** differential sensitivity state vector array */
        AmiVectorArray *sdx;
        /** Jacobian */
        DlsMat Jtmp;
        
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
        
    };
    
} // namespace amici

#endif // MODEL_H
