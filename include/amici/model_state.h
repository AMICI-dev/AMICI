#ifndef AMICI_MODEL_STATE_H
#define AMICI_MODEL_STATE_H

#include "amici/defines.h"
#include "amici/sundials_matrix_wrapper.h"
#include "amici/model_dimensions.h"

#include <vector>

namespace amici {


/**
 * @brief Exchange format to store and transfer the state of the
 * model at a specific timepoint.
 *
 * This is designed to only encompass the minimal
 * number of attributes that need to be transferred.
 */
struct ModelState {
    /**
     * Flag indicating whether a certain Heaviside function should be active or
     * not (dimension: `ne`)
     */
    std::vector<realtype> h;

    /** Total abundances for conservation laws
     (dimension: `nx_rdata - nx_solver`) */
    std::vector<realtype> total_cl;

    /** Sensitivities of total abundances for conservation laws
     (dimension: `(nx_rdata-nx_solver) x np`, row-major) */
    std::vector<realtype> stotal_cl;

    /** Unscaled parameters (dimension: `np`) */
    std::vector<realtype> unscaledParameters;

    /** Constants (dimension: `nk`) */
    std::vector<realtype> fixedParameters;

    /**
     * Indexes of parameters wrt to which sensitivities are computed
     * (dimension: nplist)
     */
    std::vector<int> plist;
};


/**
 * @brief Storage for `amici::Model` quantities computed based on
 * `amici::ModelState` for a specific timepoint.
 *
 * Serves as workspace for a model simulation to avoid repeated reallocation.
 */
struct ModelStateDerived {
    ModelStateDerived() = default;

    /**
     * @brief Constructor from model dimensions.
     * @param dim Model dimensions
     */
    explicit ModelStateDerived(ModelDimensions const& dim);

    /** Sparse Jacobian (dimension: `amici::Model::nnz`) */
    SUNMatrixWrapper J_;

    /** Sparse Backwards Jacobian (dimension: `amici::Model::nnz`) */
    SUNMatrixWrapper JB_;

    /** Sparse dxdotdw temporary storage (dimension: `ndxdotdw`) */
    SUNMatrixWrapper dxdotdw_;

    /** Sparse dwdx temporary storage (dimension: `ndwdx`) */
    SUNMatrixWrapper dwdx_;

    /** Sparse dwdp temporary storage (dimension: `ndwdp`) */
    SUNMatrixWrapper dwdp_;

    /** Dense Mass matrix (dimension: `nx_solver` x `nx_solver`) */
    SUNMatrixWrapper M_;

    /**
     * Temporary storage of `dxdotdp_full` data across functions (Python only)
     * (dimension: `nplist` x `nx_solver`, nnz: dynamic,
     * type `CSC_MAT`)
     */
    SUNMatrixWrapper dxdotdp_full;

    /**
     * Temporary storage of `dxdotdp_explicit` data across functions (Python only)
     * (dimension: `nplist` x `nx_solver`, nnz:  `ndxdotdp_explicit`,
     *  type `CSC_MAT`)
     */
    SUNMatrixWrapper dxdotdp_explicit;

    /**
     * Temporary storage of `dxdotdp_implicit` data across functions,
     * Python-only
     * (dimension: `nplist` x `nx_solver`, nnz: dynamic,
     * type `CSC_MAT`)
     */
    SUNMatrixWrapper dxdotdp_implicit;

    /**
     * Temporary storage of `dxdotdx_explicit` data across functions (Python only)
     * (dimension: `nplist` x `nx_solver`, nnz: `nxdotdotdx_explicit`,
     *  type `CSC_MAT`)
     */
    SUNMatrixWrapper dxdotdx_explicit;

    /**
     * Temporary storage of `dxdotdx_implicit` data across functions,
     * Python-only
     * (dimension: `nplist` x `nx_solver`, nnz: dynamic,
     * type `CSC_MAT`)
     */
    SUNMatrixWrapper dxdotdx_implicit;

    /**
     * Temporary storage of `dxdotdp` data across functions, Matlab only
     * (dimension: `nplist` x `nx_solver` , row-major)
     */
    AmiVectorArray dxdotdp {0, 0};

    /** Sparse observable derivative of data likelihood, only used if
     * `pythonGenerated` == `true` (dimension `nytrue`, `nJ` x `ny`, row-major)
     */
    std::vector<SUNMatrixWrapper> dJydy_;

    /** Observable derivative of data likelihood, only used if
     * `pythonGenerated` == `false` (dimension `nJ` x `ny` x `nytrue` ,
     * row-major)
     */
    std::vector<realtype> dJydy_matlab_;

    /** Observable sigma derivative of data likelihood
     * (dimension nJ x ny x nytrue, row-major)
     */
    std::vector<realtype> dJydsigma_;

    /** State derivative of data likelihood
     * (dimension `nJ` x `nx_solver`, row-major)
     */
    std::vector<realtype> dJydx_;

    /** Parameter derivative of data likelihood for current timepoint
     * (dimension: nJ x nplist, row-major)
     */
    std::vector<realtype> dJydp_;

    /** event output derivative of event likelihood
     * (dimension nJ x nz x nztrue, row-major)
     */
    std::vector<realtype> dJzdz_;

    /** event sigma derivative of event likelihood
     * (dimension nJ x nz x nztrue, row-major)
     */
    std::vector<realtype> dJzdsigma_;

    /** event output derivative of event likelihood at final timepoint
     * (dimension nJ x nz x nztrue, row-major)
     */
    std::vector<realtype> dJrzdz_;

    /** event sigma derivative of event likelihood at final timepoint
     * (dimension nJ x nz x nztrue, row-major)
     */
    std::vector<realtype> dJrzdsigma_;

    /** state derivative of event likelihood
     * (dimension `nJ` x `nx_solver`, row-major)
     */
    std::vector<realtype> dJzdx_;

    /** parameter derivative of event likelihood for current timepoint
     * (dimension: nJ x nplist x, row-major)
     */
    std::vector<realtype> dJzdp_;

    /** state derivative of event output
     * (dimension: nz x `nx_solver`, row-major)
     */
    std::vector<realtype> dzdx_;

    /** parameter derivative of event output
     * (dimension: nz x nplist, row-major)
     */
    std::vector<realtype> dzdp_;

    /** state derivative of event regularization variable
     * (dimension: `nz` x `nx_solver`, row-major)
     */
    std::vector<realtype> drzdx_;

    /** parameter derivative of event regularization variable
     * (dimension: nz x nplist, row-major)
     */
    std::vector<realtype> drzdp_;

    /** parameter derivative of observable
     * (dimension: ny x nplist, row-major)
     */
    std::vector<realtype> dydp_;

    /** state derivative of time-resolved observable
     * (dimension: `nx_solver` x `ny`, row-major)
     */
    std::vector<realtype> dydx_;

    /** temporary storage of w data across functions (dimension: nw) */
    std::vector<realtype> w_;

    /** temporary storage for flattened sx,
     * (dimension: `nx_solver` x `nplist`, row-major)
     */
    std::vector<realtype> sx_;

    /** temporary storage for `x_rdata` (dimension: `nx_rdata`) */
    std::vector<realtype> x_rdata_;

    /** temporary storage for `sx_rdata` slice (dimension: `nx_rdata`) */
    std::vector<realtype> sx_rdata_;

    /** temporary storage for time-resolved observable (dimension: ny) */
    std::vector<realtype> y_;

    /** data standard deviation for current timepoint (dimension: ny) */
    std::vector<realtype> sigmay_;

    /** temporary storage for  parameter derivative of data standard deviation,
     * (dimension: ny x nplist, row-major)
     */
    std::vector<realtype> dsigmaydp_;

    /** temporary storage for event-resolved observable (dimension: nz) */
    std::vector<realtype> z_;

    /** temporary storage for event regularization (dimension: nz) */
    std::vector<realtype> rz_;

    /** temporary storage for event standard deviation (dimension: nz) */
    std::vector<realtype> sigmaz_;

    /** temporary storage for  parameter derivative of event standard deviation,
     * (dimension: nz x nplist, row-major)
     */
    std::vector<realtype> dsigmazdp_;

    /** temporary storage for change in x after event (dimension: `nx_solver`) */
    std::vector<realtype> deltax_;

    /** temporary storage for change in sx after event
     * (dimension: `nx_solver` x `nplist`, row-major)
     */
    std::vector<realtype> deltasx_;

    /** temporary storage for change in xB after event (dimension: `nx_solver`) */
    std::vector<realtype> deltaxB_;

    /** temporary storage for change in qB after event
     * (dimension: nJ x nplist, row-major)
     */
    std::vector<realtype> deltaqB_;

    /** temporary storage of positified state variables according to
     * stateIsNonNegative (dimension: `nx_solver`) */
    AmiVector x_pos_tmp_ {0};
};


} // namespace amici

#endif // AMICI_MODEL_STATE_H
