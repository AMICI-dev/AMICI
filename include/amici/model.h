#ifndef AMICI_MODEL_H
#define AMICI_MODEL_H

#include "amici/exception.h"
#include "amici/defines.h"
#include "amici/vector.h"
#include "amici/symbolic_functions.h"

#include <sundials/sundials_direct.h> // DlsMat
#include <sundials/sundials_sparse.h> // SlsMat

#include <numeric>
#include <vector>
#include <memory>

namespace amici {

class ReturnData;
class ExpData;
class Model;
class Solver;

} // namespace amici


// for serialization friend in amici::Model
namespace boost { namespace serialization {
template <class Archive>
void serialize(Archive &ar, amici::Model &u, const unsigned int version);
}}


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
        : nx(0), nxtrue(0), ny(0), nytrue(0), nz(0), nztrue(0),
        ne(0), nw(0), ndwdx(0), ndwdp(0), nnz(0), nJ(0), ubw(0), lbw(0),
        o2mode(SecondOrderMode::none) {}
        
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
              amici::SecondOrderMode o2mode, const std::vector<amici::realtype> &p,
              std::vector<amici::realtype> k, const std::vector<int> &plist,
              std::vector<amici::realtype> idlist, std::vector<int> z2event);
        
        /** Copy constructor
         * @param other object to copy from
         * @return
         */
        Model(Model const& other);

        /** destructor */
        virtual ~Model();

        /** Copy assignment is disabled until const members are removed
         * @param other object to copy from
         * @return
         */
        Model& operator=(Model const &other)=delete;

        /**
         * @brief Clone this instance
         * @return The clone
         */
        virtual Model* clone() const = 0;

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
        
        /** Initial states
         * @param x pointer to state variables
         */
        void fx0(AmiVector *x);
        
        /** Sets only those initial states that are specified via fixedParmeters
         * @param x pointer to state variables
         */
        void fx0_fixedParameters(AmiVector *x);
        
        /** Initial value for time derivative of states (only necessary for DAEs)
         * @param x0 Vector with the initial states
         * @param dx0 Vector to which the initial derivative states will be
         * written (only DAE)
         **/
        virtual void fdx0(AmiVector *x0, AmiVector *dx0);

        /** Initial value for initial state sensitivities
          * @param sx pointer to state sensitivity variables
          * @param x pointer to state variables
          **/

        void fsx0(AmiVectorArray *sx, const AmiVector *x);
        
        /** Sensitivity of derivative initial states sensitivities sdx0 (only
         *  necessary for DAEs)
         **/
        virtual void fsdx0();
        
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
        
        void fsigmay(const int it, ReturnData *rdata, const ExpData *edata);
        
        void fdsigmaydp(const int it, ReturnData *rdata, const ExpData *edata);
        
        void fsigmaz(const realtype t, const int ie, const int *nroots, ReturnData *rdata,
                     const ExpData *edata);
        
        void fdsigmazdp(const realtype t, const int ie, const int *nroots, ReturnData *rdata, const ExpData *edata);
        
        void fJy(const int it, ReturnData *rdata, const ExpData *edata);
        
        void fJz(const int nroots, ReturnData *rdata, const ExpData *edata);
        
        void fJrz(const int nroots, ReturnData *rdata, const ExpData *edata);

        void fdJydy(const int it, const ReturnData *rdata, const ExpData *edata);
        
        void fdJydsigma(const int it, const ReturnData *rdata, const ExpData *edata);
        
        void fdJzdz(const int nroots, const ReturnData *rdata, const ExpData *edata);
        
        void fdJzdsigma(const int nroots, const ReturnData *rdata, const ExpData *edata);
        
        void fdJrzdz(const int nroots, const ReturnData *rdata, const ExpData *edata);
        
        void fdJrzdsigma(const int nroots,const ReturnData *rdata, const ExpData *edata);
                
        // Generic implementations
        void fsy(const int it, ReturnData *rdata);
        
        void fsz_tf(const int *nroots, const int ie, ReturnData *rdata);
        
        void fsJy(const int it, const std::vector<realtype>& dJydx, ReturnData *rdata);
        
        void fdJydp(const int it, const ExpData *edata,
                   ReturnData *rdata);
        
        void fdJydx(std::vector<realtype> *dJydx, const int it, const ExpData *edata, const ReturnData *rdata);
        
        void fsJz(const int nroots, const std::vector<realtype>& dJzdx, AmiVectorArray *sx, ReturnData *rdata);
        
        void fdJzdp(const int nroots, realtype t, const ExpData *edata, const ReturnData *rdata);
        
        void fdJzdx(std::vector<realtype> *dJzdx, const int nroots, realtype t, const ExpData *edata, const ReturnData *rdata);
        
        void initialize(AmiVector *x, AmiVector *dx);
        
        void initializeStates(AmiVector *x);
        
        void initHeaviside(AmiVector *x, AmiVector *dx);
        
        /** number of paramaeters wrt to which sensitivities are computed
         * @return length of sensitivity index vector
         */
        int nplist() const;

        /** total number of model parameters
         * @return length of parameter vector
         */
        int np() const;

        /** number of constants
         * @return length of constant vector
         */
        int nk() const;

        /** constants
         * @return pointer to constants array
         */
        const double *k() const;

        /**
         * @brief Get nmaxevent
         * @return nmaxevent
         */
        int nMaxEvent() const;

        /**
         * @brief setNMaxEvent
         * @param nmaxevent
         */
        void setNMaxEvent(int nmaxevent);

        /**
         * @brief Get number of timepoints
         * @return number of timepoints
         */
        int nt() const;

        /**
         * @brief getParameterScale
         * @return
         */
        std::vector<ParameterScaling> const& getParameterScale() const;

        /**
         * @brief setParameterScale
         * @param pscale
         */
        void setParameterScale(ParameterScaling pscale);

        /**
         * @brief setParameterScale
         * @param pscale
         */
        void setParameterScale(const std::vector<ParameterScaling>& pscale);

        /**
         * @brief getParameters
         * @return The user-set parameters (see also getUnscaledParameters)
         */
        std::vector<realtype> const& getParameters() const;

        /**
         * @brief setParameters
         * @param p
         */
        void setParameters(std::vector<realtype> const&p);

        /**
         * @brief getUnscaledParameters
         * @return The unscaled parameters
         */
        std::vector<realtype> const& getUnscaledParameters() const;

        /**
         * @brief getFixedParameters
         * @return
         */
        std::vector<realtype> const& getFixedParameters() const;

        /**
         * @brief setFixedParameters
         * @param k
         */
        void setFixedParameters(std::vector<realtype> const&k);

        /**
         * @brief getTimepoints
         * @return
         */
        std::vector<realtype> const& getTimepoints() const;

        /**
         * @brief setTimepoints
         * @param ts
         */
        void setTimepoints(std::vector<realtype> const& ts);

        /**
         * @brief Get timepoint for given index
         * @param idx
         * @return
         */
        double t(int idx) const;

        /**
         * @brief getParameterList
         * @return
         */
        std::vector<int> const& getParameterList() const;

        /**
         * @brief setParameterList
         * @param plist
         */
        void setParameterList(std::vector<int> const& plist);

        /**
         * @brief getInitialStates
         * @return
         */
        std::vector<realtype> const& getInitialStates() const;

        /**
         * @brief setInitialStates
         * @param x0
         */
        void setInitialStates(std::vector<realtype> const& x0);

        /**
         * @brief getInitialStateSensitivities
         * @return
         */
        std::vector<realtype> const& getInitialStateSensitivities() const;

        /**
         * @brief setInitialStateSensitivities
         * @param sx0
         */
        void setInitialStateSensitivities(std::vector<realtype> const& sx0);

        /** initial timepoint
         *  @return timepoint
         */
        double t0() const;

        /**
         * @brief setT0
         * @param t0
         */
        void setT0(double t0);

        /** entry in parameter list
          * @param pos index
          * @return entry
          */
        int plist(int pos) const;

        /**
         * @brief unscaleParameters
         * @param bufferUnscaled
         */
        void unscaleParameters(double *bufferUnscaled) const;

        /**
         * @brief Require computation of sensitivities for all parameters p [0..np[
         * in natural order.
         */
        void requireSensitivitiesForAllParameters();

        void fw(const realtype t, const N_Vector x);

        void fdwdp(const realtype t, const N_Vector x);
        
        void fdwdx(const realtype t, const N_Vector x);
        
        /** residual function
         * @param it time index
         * @param rdata ReturnData instance to which result will be written
         * @param edata ExpData instance containing observable data
         */
        void fres(const int it, ReturnData *rdata, const ExpData *edata);
        
        /** chi-squared function
         * @param it time index
         * @param rdata ReturnData instance to which result will be written
         */
        void fchi2(const int it, ReturnData *rdata);
        
        /** residual sensitivity function
         * @param it time index
         * @param rdata ReturnData instance to which result will be written
         * @param edata ExpData instance containing observable data
         */
        void fsres(const int it, ReturnData *rdata, const ExpData *edata);
        
        /** fisher information matrix function
         * @param it time index
         * @param rdata ReturnData instance to which result will be written
         */
        void fFIM(const int it, ReturnData *rdata);

        
        /**
         * updateHeaviside updates the heaviside variables h on event occurences
         *
         * @param rootsfound provides the direction of the zero-crossing, so adding
         it will give the right update to the heaviside variables (zero if no root
         was found)
         */
        void updateHeaviside(const std::vector<int>& rootsfound);
        
        /**
         * updateHeavisideB updates the heaviside variables h on event occurences
         in the backward problem
         *
         * @param rootsfound provides the direction of the zero-crossing, so adding
         it will give the right update to the heaviside variables (zero if no root
         was found)
         */
        void updateHeavisideB(const int *rootsfound);
        
        /**
         * @brief Serialize Model (see boost::serialization::serialize)
         * @param ar Archive to serialize to
         * @param r Data to serialize
         * @param version Version number
         */
        template <class Archive>
        friend void boost::serialization::serialize(Archive &ar, Model &r, const unsigned int version);

        /**
         * @brief Check equality of data members
         * @param a
         * @param b
         * @return
         */
        friend bool operator ==(const Model &a, const Model &b);
        
        realtype gett(const int it, const ReturnData *rdata) const;

        /**
         * @brief Check if the given array has only finite elements.
         * If not try to give hints by which other fields this could be caused.
         * @param N
         * @param array
         * @param fun
         * @return AMICI_RECOVERABLE_ERROR if a NaN/Inf value was found, AMICI_SUCCESS otherwise
         */
        int checkFinite(const int N,const realtype *array, const char* fun) const;


        /**
         * @brief Reports whether the model has parameter names set.
         * @return
         */
        virtual bool hasParameterNames() const { return np() && !getParameterNames().empty(); }

        /**
         * @brief Get names of the model parameters
         * @return the names
         */
        virtual std::vector<std::string> getParameterNames() const { return std::vector<std::string>(); }

        /**
         * @brief Reports whether the model has state names set.
         * @return
         */
        virtual bool hasStateNames() const { return nx && !getStateNames().empty(); }

        /**
         * @brief Get names of the model states
         * @return the names
         */
        virtual std::vector<std::string> getStateNames() const { return std::vector<std::string>(); }

        /**
         * @brief Reports whether the model has fixed parameter names set.
         * @return
         */
        virtual bool hasFixedParameterNames() const { return nk() && !getFixedParameterNames().empty(); }

        /**
         * @brief Get names of the fixed model parameters
         * @return the names
         */
        virtual std::vector<std::string> getFixedParameterNames() const { return std::vector<std::string>(); }

        /**
         * @brief Reports whether the model has observable names set.
         * @return
         */
        virtual bool hasObservableNames() const { return ny && !getObservableNames().empty(); }

        /**
         * @brief Get names of the observables
         * @return the names
         */
        virtual std::vector<std::string> getObservableNames() const { return std::vector<std::string>(); }
        
        /**
         * @brief Reports whether the model has parameter ids set.
         * @return
         */
        virtual bool hasParameterIds() const { return np() && !getParameterIds().empty(); }
        
        /**
         * @brief Get ids of the model parameters
         * @return the ids
         */
        virtual std::vector<std::string> getParameterIds() const { return std::vector<std::string>(); }
        
        /**
         * @brief Reports whether the model has state ids set.
         * @return
         */
        virtual bool hasStateIds() const { return nx && !getStateIds().empty(); }
        
        /**
         * @brief Get ids of the model states
         * @return the ids
         */
        virtual std::vector<std::string> getStateIds() const { return std::vector<std::string>(); }
        
        /**
         * @brief Reports whether the model has fixed parameter ids set.
         * @return
         */
        virtual bool hasFixedParameterIds() const { return nk() && !getFixedParameterIds().empty(); }
        
        /**
         * @brief Get ids of the fixed model parameters
         * @return the ids
         */
        virtual std::vector<std::string> getFixedParameterIds() const { return std::vector<std::string>(); }
        
        /**
         * @brief Reports whether the model has observable ids set.
         * @return
         */
        virtual bool hasObservableIds() const { return ny && !getObservableIds().empty(); }
        
        /**
         * @brief Get ids of the observables
         * @return the ids
         */
        virtual std::vector<std::string> getObservableIds() const { return std::vector<std::string>(); }
        
        /**
         * @brief sets the mode how sensitivities are computed in the steadystate simulation
         * @param mode steadyStateSensitivityMode
         */
        void setSteadyStateSensitivityMode (const SteadyStateSensitivityMode mode) {
            steadyStateSensitivityMode = mode;
        }
        
        /**
         * @brief gets the mode how sensitivities are computed in the steadystate simulation
         * @return flag value
         */
        SteadyStateSensitivityMode getSteadyStateSensitivityMode () const {
            return steadyStateSensitivityMode;
        }

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
        const SecondOrderMode o2mode;
        /** index indicating to which event an event output belongs */
        const std::vector<int> z2event;
        /** flag array for DAE equations */
        const std::vector<realtype> idlist;

        /** data standard deviation for current timepoint (dimension: ny) */
        std::vector<realtype> sigmay;
        /** parameter derivative of data standard deviation for current timepoint (dimension: nplist x ny, row-major) */
        std::vector<realtype> dsigmaydp;
        /** event standard deviation for current timepoint (dimension: nz) */
        std::vector<realtype> sigmaz;
        /** parameter derivative of event standard deviation for current timepoint (dimension: nplist x nz, row-major) */
        std::vector<realtype> dsigmazdp;
        /** parameter derivative of data likelihood for current timepoint (dimension: nplist x nJ, row-major) */
        std::vector<realtype> dJydp;
        /** parameter derivative of event likelihood for current timepoint (dimension: nplist x nJ, row-major) */
        std::vector<realtype> dJzdp;

        /** change in x at current timepoint (dimension: nx) */
        std::vector<realtype> deltax;
        /** change in sx at current timepoint (dimension: nplist x nx, row-major) */
        std::vector<realtype> deltasx;
        /** change in xB at current timepoint (dimension: nJ x nxtrue, row-major) */
        std::vector<realtype> deltaxB;
        /** change in qB at current timepoint (dimension: nJ x nplist, row-major) */
        std::vector<realtype> deltaqB;

        /** tempory storage of dxdotdp data across functions (dimension: nplist x nx, row-major) */
        std::vector<realtype> dxdotdp;


    protected:

        /**
         * @brief Set the nplist-dependent vectors to their proper sizes
         */
        void initializeVectors();

        
        /** model specific implementation of fx0
         * @param x0 initial state
         * @param t initial time
         * @param p parameter vector
         * @param k constant vector
         **/
        virtual void fx0(realtype *x0, const realtype t, const realtype *p, const realtype *k) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fx0
         * @param x0 initial state
         * @param t initial time
         * @param p parameter vector
         * @param k constant vector
         **/
        virtual void fx0_fixedParameters(realtype *x0, const realtype t, const realtype *p, const realtype *k) {
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
        
        /** model specific implementation of fdrzdp
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
        
        /** model specific implementation of fdeltaqB
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
        virtual void fsigmay(realtype *sigmay, const realtype t, const realtype *p, const realtype *k) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fsigmay
         * @param dsigmaydp partial derivative of standard deviation of measurements
         * @param t current time
         * @param p parameter vector
         * @param k constant vector
         * @param ip sensitivity index
         **/
        virtual void fdsigmaydp(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const int ip) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fsigmaz
         * @param sigmaz standard deviation of event measurements
         * @param t current time
         * @param p parameter vector
         * @param k constant vector
         **/
        virtual void fsigmaz(realtype *sigmaz, const realtype t, const realtype *p, const realtype *k) {
            throw AmiException("Requested functionality is not supported as (%s) is not implemented for this model!",__func__);
        }
        
        /** model specific implementation of fsigmaz
         * @param dsigmazdp partial derivative of standard deviation of event measurements
         * @param t current time
         * @param p parameter vector
         * @param k constant vector
         * @param ip sensitivity index
         **/
        virtual void fdsigmazdp(realtype *dsigmazdp, const realtype t, const realtype *p, const realtype *k, const int ip) {
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
                        const realtype *k, const realtype *h) {}
        
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
                           const realtype *k, const realtype *h, const realtype *w) {}
        
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
                           const realtype *k, const realtype *h, const realtype *w) {}
        
        void getmy(const int it, const ExpData *edata);
        
        void getmz(const int nroots, const ExpData *edata);
        
        const realtype *gety(const int it, const ReturnData *rdata) const;
        
        const realtype *getx(const int it, const ReturnData *rdata) const;
        
        const realtype *getsx(const int it, const ReturnData *rdata) const;
        
        const realtype *getz(const int nroots, const ReturnData *rdata) const;
        
        const realtype *getrz(const int nroots, const ReturnData *rdata) const;
        
        const realtype *getsz(const int nroots, const int nplist, const ReturnData *rdata) const;
        
        const realtype *getsrz(const int nroots, const int nplist, const ReturnData *rdata) const;



        /** Sparse Jacobian (dimension: nnz)*/
        SlsMat J = nullptr;
        
        /** current observable (dimension: nytrue) */
        std::vector<realtype> my;
        /** current event measurement (dimension: nztrue) */
        std::vector<realtype> mz;
        /** observable derivative of data likelihood (dimension nJ x nytrue x ny, ordering = ?) */
        std::vector<realtype> dJydy;
        /** observable sigma derivative of data likelihood (dimension nJ x nytrue x ny, ordering = ?) */
        std::vector<realtype> dJydsigma;

        /** event ouput derivative of event likelihood (dimension nJ x nztrue x nz, ordering = ?) */
        std::vector<realtype> dJzdz;
        /** event sigma derivative of event likelihood (dimension nJ x nztrue x nz, ordering = ?) */
        std::vector<realtype> dJzdsigma;
        /** event ouput derivative of event likelihood at final timepoint (dimension nJ x nztrue x nz, ordering = ?) */
        std::vector<realtype> dJrzdz;
        /** event sigma derivative of event likelihood at final timepoint (dimension nJ x nztrue x nz, ordering = ?) */
        std::vector<realtype> dJrzdsigma;
        /** state derivative of event output (dimension: nz * nx, ordering = ?) */
        std::vector<realtype> dzdx;
        /** parameter derivative of event output (dimension: nz * nplist, ordering = ?) */
        std::vector<realtype> dzdp;
        /** state derivative of event timepoint (dimension: nz * nx, ordering = ?) */
        std::vector<realtype> drzdx;
        /** parameter derivative of event timepoint (dimension: nz * nplist, ordering = ?) */
        std::vector<realtype> drzdp;
        /** parameter derivative of observable (dimension: nplist * ny, row-major) */
        std::vector<realtype> dydp;

        /** state derivative of observable (dimension: ny * nx, ordering = ?) */
        std::vector<realtype> dydx;
        /** tempory storage of w data across functions (dimension: nw) */
        std::vector<realtype> w;
        /** tempory storage of sparse dwdx data across functions (dimension: ndwdx) */
        std::vector<realtype> dwdx;
        /** tempory storage of sparse dwdp data across functions (dimension: ndwdp) */
        std::vector<realtype> dwdp;
        /** tempory storage of M data across functions (dimension: nx) */
        std::vector<realtype> M;
        /** tempory storage of stau data across functions (dimension: nplist) */
        std::vector<realtype> stau;
        
        /** flag indicating whether a certain heaviside function should be active or
         not (dimension: ne) */
        std::vector<realtype> h;

        /** unscaled parameters (dimension: np) */
        std::vector<realtype> unscaledParameters;

        /** orignal user-provided, possibly scaled parameter array (size np) */
        std::vector<realtype>originalParameters;

        /** constants (dimension: nk) */
        std::vector<realtype> fixedParameters;

        /** indexes of parameters wrt to which sensitivities are computed (dimension nplist) */
        std::vector<int> plist_;

        /** state initialisation (size nx) */
        std::vector<double>x0data;

        /** sensitivity initialisation (size nx * nplist, ordering = ?) */
        std::vector<realtype>sx0data;

        /** timepoints (size nt) */
        std::vector<realtype>ts;

        /** maximal number of events to track */
        int nmaxevent = 10;

        /** parameter transformation of `originalParameters` (dimension np) */
        std::vector<ParameterScaling> pscale;

        /** starting time */
        double tstart = 0.0;
        
        /** flag indicating whether steadystate sensivities are to be computed
         via FSA when steadyStateSimulation is used */
        SteadyStateSensitivityMode steadyStateSensitivityMode = SteadyStateSensitivityMode::newtonOnly;

    };

    bool operator ==(const Model &a, const Model &b);
    
} // namespace amici

#endif // AMICI_MODEL_H
