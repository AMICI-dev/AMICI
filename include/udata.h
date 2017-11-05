#ifndef AMICI_UDATA_H
#define AMICI_UDATA_H

#include "include/amici_defines.h"
#include <cmath>
#include <vector>

namespace amici {
class UserData;
}

namespace boost { namespace serialization {
template <class Archive>
void serialize(Archive &ar, amici::UserData &u, const unsigned int version);
}}

namespace amici {

/** @brief struct that stores all user provided data
 * NOTE: multidimensional arrays are expected to be stored in column-major order
 * (FORTRAN-style)
 */
class UserData {

  public:
    /**
     * @brief UserData
     * @param np total number of model parameters
     * @param nk number of fixed parameters
     * @param nx number of states
     */
    UserData(int np, int nk, int nx);

    /**
     * @brief Default constructor for testing and serialization
     */
    UserData();

    /**
     * @brief Copy constructor
     * @param other object to copy from
     */
    UserData (const UserData &other);

    /**
     * @brief Copy assignment is disabled until const members are removed
     * @param other object to copy from
     * @return
     */
    UserData& operator=(UserData const &other)=delete;

    void unscaleParameters(double *bufferUnscaled) const;

    /**
     * @brief setTimepoints
     * @param timepoints
     * @param numTimepoints
     */
    void setTimepoints(const double *timepoints, int numTimepoints);

    /**
     * @brief setParameters
     * @param parameters
     */
    void setParameters(const double *parameters);

    /**
     * @brief setConstants
     * @param constants
     */
    void setConstants(const double *constants);

    /**
     * @brief setPlist set parameter selection and ordering.
     * accepts array of doubles to deal with user input from matlab.
     * @param plist
     * @param nplist length of plist
     */
    void setPlist(const double *plist, int nplist);
    /**
     * @brief setPlist set parameter selection and ordering.
     * accepts array of ints.
     * @param plist
     * @param nplist length of plist
     */
    void setPlist(const int *plist, int nplist);

    /**
     * @brief Require computation of sensitiivities for all parameters p [0..np[
     * in natural order.
     */
    void requireSensitivitiesForAllParameters();

    /**
     * @brief setPbar. Must not be called before setPlist
     * @param parameterScaling
     */
    void setPbar(const double *parameterScaling);
    
    /**
     * @brief setXbar.
     * @param stateScaling
     */
    void setXbar(const double *stateScaling);

    /**
     * @brief setStateInitialization
     * @param stateInitialization
     */
    void setStateInitialization(const double *stateInitialization);

    /**
     * @brief setSensitivityInitialization
     * @param sensitivityInitialization
     */
    void setSensitivityInitialization(const double *sensitivityInitialization);

    ~UserData();

    /* Options */


    /** number of parameters in plist 
     *  @return length of plist
     */
    const int nplist() const{
        return p_index.size();
    };
    
    const double tstart() const{
        return t0;
    }

    /** number of timepoints */
    const int nt() const{
        return ts.size();
    }

    /** function to print the contents of the UserData object */
    void print() const;

    /** total number of model parameters */
    const int np() const{
        return par.size();
    };
    /** number of fixed parameters */
    const int nk() const{
        return konst.size();
    };
    /** number of states */
    const int nx() const{
        return sizex;
    };
    
    /** max number of events */
    const int nme() const{
        return nmaxevent;
    };
    
    /** max number of events */
    const int sensmeth() const{
        return sensi_meth;
    };
    
    const double *p() const{
        return par.data();
    };
    
    const double *k() const{
        return konst.data();
    };
    
    const int plist(int pos) const{
        return p_index.at(pos);
    };
    
    const std::vector<double> getInitialStates() const {
        return x0data;
    }

    /**
     * @brief Serialize UserData (see boost::serialization::serialize)
     * @param ar Archive to serialize to
     * @param r Data to serialize
     * @param version Version number
     */
    template <class Archive>
    friend void boost::serialization::serialize(Archive &ar, UserData &r, const unsigned int version);

private:
    /** maximal number of events to track */
    int nmaxevent = 10;
    
    /** positivity flag (size nx) */
    std::vector<int>qpositivex;
    
    /** parameter selection and reordering (size nplist) */
    std::vector<int>p_index;
    
    
    const int sizex;
    /** parameter array (size np) */
    std::vector<double>par;
    
    /** constants array (size nk) */
    std::vector<double>konst;
    
    /** starting time */
    double t0 = 0.0;
    
    /** parameter transformation of p */
    AMICI_parameter_scaling pscale = AMICI_SCALING_NONE;
    
    /** timepoints (size nt) */
    std::vector<double>ts;
    
    /** scaling of parameters (size nplist) */
    std::vector<double>pbar;
    
    /** scaling of states (size nx)
     * NOTE: currently not used */
    std::vector<double>xbar;
    
    /** flag indicating whether sensitivities are supposed to be computed */
    AMICI_sensi_order sensi = AMICI_SENSI_ORDER_NONE;
    
    /** absolute tolerances for integration */
    double atol = 1e-16;
    
    /** relative tolerances for integration */
    double rtol = 1e-8;
    
    /** maximum number of allowed integration steps */
    int maxsteps = 0;
    
    /** maximum number of allowed Newton steps for steady state computation */
    int newton_maxsteps = 0;
    
    /** maximum number of allowed linear steps per Newton step for steady state
     * computation */
    int newton_maxlinsteps = 0;
    
    /** Preequilibration of model via Newton solver? */
    int newton_preeq = false;
    
    /** Which preconditioner is to be used in the case of iterative linear
     * Newton solvers */
    int newton_precon = 1;
    
    /** internal sensitivity method flag used to select the sensitivity solution
     * method. Only applies for Forward Sensitivities. */
    InternalSensitivityMethod ism = SIMULTANEOUS;
    
    /** method for sensitivity computation */
    AMICI_sensi_meth sensi_meth = AMICI_SENSI_FSA;
    
    /** linear solver specification */
    LinearSolver linsol = AMICI_KLU;
    
    /** interpolation type for the forward problem solution which
     * is then used for the backwards problem.
     */
    InterpolationType interpType = HERMITE;
    
    /** specifies the linear multistep method.
     */
    LinearMultistepMethod lmm = BDF;
    
    /**
     * specifies the type of nonlinear solver iteration
     */
    NonlinearSolverIteration iter = NEWTON;
    
    /** flag controlling stability limit detection */
    booleantype stldet = true;
    
    /** state initialisation (size nx) */
    std::vector<double>x0data;
    
    /** sensitivity initialisation (size nx * nplist) */
    std::vector<double>sx0data;
    
    /** state ordering */
    StateOrdering ordering = AMD;
};

} // namespace amici

#endif /* _MY_UDATA */
