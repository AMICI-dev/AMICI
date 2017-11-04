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
    class Model_ODE {
    public:
        /** default constructor */
        Model_ODE() : Model() {}
        
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
        Model_ODE(const int np, const int nx, const int nxtrue, const int nk,
              const int ny, const int nytrue, const int nz, const int nztrue,
              const int ne, const int nJ, const int nw, const int ndwdx,
              const int ndwdp, const int nnz, const int ubw, const int lbw,
              const AMICI_o2mode o2mode)
        : Model(np,nx,nxtru,nk,ny,nytrue,nz,nztrue,ne,nJ,nw,ndwdx,ndwdp,nnz,ubw,lbw,o2mode){}
        
    };
    
} // namespace amici

#endif // MODEL_H
