#ifndef AMICI_MODEL_DIMENSIONS_H
#define AMICI_MODEL_DIMENSIONS_H

#include <gsl/gsl-lite.hpp>
#include <vector>

namespace amici {

/**
 * @brief Container for model dimensions.
 *
 * Holds number of state variables, observables, etc.
 */
struct ModelDimensions {
    /**
     * @brief Validate dimensions.
     */
    void validate() const {
        Expects(nxtrue_rdata >= 0);
        Expects(nxtrue_rdata <= nx_rdata);
        Expects(nxtrue_solver >= 0);
        Expects(nx_solver <= nx_rdata);
        Expects(nxtrue_solver <= nx_solver);
        Expects(nx_solver_reinit >= 0);
        Expects(nx_solver_reinit <= nx_solver);
        Expects(np >= 0);
        Expects(nk >= 0);
        Expects(nytrue <= ny);
        Expects(nytrue >= 0);
        Expects(nztrue >= 0);
        Expects(nztrue <= nz);
        Expects(ne >= 0);
        Expects(ne_solver >= 0);
        Expects(ne >= ne_solver);
        Expects(nspl >= 0);
        Expects(nw >= 0);
        Expects(ndwdx >= 0);
        Expects(ndwdx <= nw * nx_solver);
        Expects(ndwdp >= 0);
        Expects(ndwdp <= nw * np);
        Expects(ndwdw >= 0);
        Expects(ndwdw <= nw * nw);
        Expects(ndxdotdw >= 0);
        Expects(ndxrdatadxsolver >= 0);
        Expects(ndxrdatadxsolver <= nx_rdata * nx_solver);
        Expects(ndxrdatadtcl >= 0);
        Expects(ndxrdatadtcl <= nx_rdata * (nx_rdata - nx_solver));
        Expects(ndtotal_cldx_rdata >= 0);
        Expects(ndtotal_cldx_rdata <= (nx_rdata - nx_solver) * nx_rdata);
        Expects(nnz >= 0);
        Expects(nJ >= 0);
        Expects(ubw >= 0);
        Expects(lbw >= 0);
        Expects(ndxdotdp_explicit >= 0);
        Expects(ndxdotdx_explicit >= 0);
        Expects(w_recursion_depth >= 0);
    }
    /** Number of state variables */
    int nx_rdata{0};

    /** Number of state variables in the unaugmented system */
    int nxtrue_rdata{0};

    /** Number of state variables with conservation laws applied */
    int nx_solver{0};

    /**
     * Number of state variables in the unaugmented system with conservation
     * laws applied
     */
    int nxtrue_solver{0};

    /** Number of solver state variables subject to reinitialization */
    int nx_solver_reinit{0};

    /** Number of parameters */
    int np{0};

    /** Number of constants */
    int nk{0};

    /** Number of observables */
    int ny{0};

    /** Number of observables in the unaugmented system */
    int nytrue{0};

    /** Number of event outputs */
    int nz{0};

    /** Number of event outputs in the unaugmented system */
    int nztrue{0};

    /** Number of events */
    int ne{0};

    /** Number of events that require root-finding */
    int ne_solver{0};

    /** Number of spline functions in the model */
    int nspl{0};

    /** Number of common expressions */
    int nw{0};

    /**
     * Number of nonzero elements in the `x` derivative of the
     * repeating elements
     */
    int ndwdx{0};

    /**
     * Number of nonzero elements in the `p` derivative of the
     * repeating elements
     */
    int ndwdp{0};

    /**
     * Number of nonzero elements in the `w` derivative of the
     * repeating elements
     */
    int ndwdw{0};

    /** Number of nonzero elements in the \f$ w \f$ derivative of \f$ xdot \f$
     */
    int ndxdotdw{0};

    /**
     * Number of nonzero elements in the \f$ y \f$ derivative of
     * \f$ dJy \f$ (dimension `nytrue`)
     */
    std::vector<int> ndJydy;

    /** Number of nonzero elements in the \f$ x \f$ derivative of \f$ x_rdata
     * \f$ */
    int ndxrdatadxsolver{0};

    /** Number of nonzero elements in the \f$ tcl\f$ derivative of \f$ x_rdata
     * \f$ */
    int ndxrdatadtcl{0};

    /** Number of nonzero elements in the \f$ x_rdata\f$ derivative of
     *  \f$ total_cl \f$ */
    int ndtotal_cldx_rdata{0};

    /** Number of nonzero entries in Jacobian */
    int nnz{0};

    /** Dimension of the augmented objective function for 2nd order ASA */
    int nJ{0};

    /** Upper bandwidth of the Jacobian */
    int ubw{0};

    /** Lower bandwidth of the Jacobian */
    int lbw{0};

    /** Number of nonzero elements in `dxdotdp_explicit` */
    int ndxdotdp_explicit = 0;

    /** Number of nonzero elements in `dxdotdx_explicit` */
    int ndxdotdx_explicit = 0;

    /** Recursion depth of fw */
    int w_recursion_depth = 0;
};

} // namespace amici

#endif // AMICI_MODEL_DIMENSIONS_H
