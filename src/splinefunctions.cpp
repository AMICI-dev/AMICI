#include "amici/splinefunctions.h"
#include "amici/defines.h"
#include "amici/vector.h"
#include "amici/amici.h"
#include "amici/exception.h"
#include <vector>


namespace amici {

static realtype evaluatePolynomial(realtype x, realtype *coeff) {
    /* Use Horner's method for nuermical efficiency:
     * spline(t) = a * t**3 + b * t**2 + c * t + d
     *           = d + t * (c + t * (b + t * a))
     * with coeff[0, 1, 2, 3] = [d, c, b, a]
     * */

    return coeff[0] + x * (coeff[1] + x * (coeff[2] + x * coeff[3]));
}

AbstractSpline::AbstractSpline(std::vector<realtype> nodes,
                               std::vector<realtype> node_values,
                               bool equidistant_spacing,
                               bool logarithmic_paraterization)
    : nodes_(nodes), node_values_(node_values),
    equidistant_spacing_(equidistant_spacing),
    logarithmic_paraterization_(logarithmic_paraterization) {

    /* we want to set the number of nodes */
    n_nodes_ = node_values.size();

    /* In case we have equidistant spacing, compute node locations */
    if (equidistant_spacing_) {
        if (nodes.size() != 2)
            throw AmiException("Splines with equidistant spacing need a nodes "
                               "vector with two elements (first/last node).");
        realtype node_start = nodes[0];
        realtype node_step = (nodes[1] - nodes[0]) / (n_nodes_ - 1);
        nodes_.resize(n_nodes_);
        nodes_[n_nodes_ - 1] = nodes[1];
        for (int i_node = 0; i_node < n_nodes_ - 1; i_node++)
            nodes_[i_node] = node_start + i_node * node_step;
    }

    if (logarithmic_paraterization_)
        for (int iNode = 0; iNode < n_nodes_; iNode++)
            node_values_[iNode] = std::log(node_values_);
}

realtype AbstractSpline::getFinalValue(){
    return finalValue_;
}

void AbstractSpline::setFinalValue(realtype finalValue){
    finalValue_ = finalValue;
}

realtype AbstractSpline::getFinalSensitivity(const int ip){
    return finalSensitivity_[ip];
}

void AbstractSpline::setFinalSensitivity(std::vector<realtype> const finalSensitivity){
    finalSensitivity_ = finalSensitivity;
}

bool AbstractSpline::get_equidistant_spacing() {
    return equidistant_spacing_;
}

void AbstractSpline::set_equidistant_spacing(bool equidistant_spacing) {
    equidistant_spacing_ = equidistant_spacing;
}

bool get_logarithmic_paraterization() {
    return logarithmic_paraterization_;
}

void AbstractSpline::set_logarithmic_paraterization(bool logarithmic_paraterization) {
    logarithmic_paraterization_ = logarithmic_paraterization;
}


HermiteSpline::HermiteSpline(std::vector<realtype> nodes,
                             std::vector<realtype> node_values,
                             std::vector<realtype> node_values_derivative,
                             SplineBoundaryCondition firstNodeBC,
                             SplineBoundaryCondition lastNodeBC,
                             SplineExtrapolation firstNodeExtrapol,
                             SplineExtrapolation lastNodeExtrapol,
                             bool node_derivative_by_FD,
                             bool equidistant_spacing,
                             bool logarithmic_paraterization)
    : AbstractSpline(nodes, node_values, equidistant_spacing,
                     logarithmic_paraterization),
    node_values_derivative_(node_values_derivative),
    firstNodeBC_(firstNodeBC), lastNodeBC_(lastNodeBC),
    firstNodeEP_(firstNodeDerivative),
    lastNodeEP_(lastNodeDerivative),
    node_derivative_by_FD_(node_derivative_by_FD) {

    /* We may have to compute the derivatives at the nodes */
    handleInnerDerviatives();
    /* First and last node need to be handeled separately */
    handleBoundaryConditions();
}

void HermiteSpline::handleInnerDerviatives() {
    /* If values of the derivative at the nodes are to be computed by finite
     * differences, we have to fill up node_values_derivative_ */
    if (node_derivative_by_FD_) {
        node_values_derivative_.resize(n_nodes(), 0.0);
        for (int i_node = 1; i_node < n_nodes() - 1; i_node++)
            node_values_derivative_[i_node] =
                (node_values_[i_node + 1] - node_values_[i_node - 1]) /
                (nodes_[i_node + 1] - nodes_[i_node - 1]);
    }
}

void HermiteSpline::handleBoundaryConditions() {
    int last = n_nodes() - 1;

    /* We have to take special care of the first node */
    switch (firstNodeBC_) {
        case SplineBoundaryCondition::given:
            if (node_derivative_by_FD_)
                /* 1-sided FD */
                node_values_derivative_[0] =
                    (node_values_[1] - node_values_[0]) /
                    (nodes_[1] - nodes_[0]);
            break;

        case SplineBoundaryCondition::zeroDerivative:
            node_values_derivative_[0] = 0;
            break;

        case SplineBoundaryCondition::natural:
            node_values_derivative_[0] = -0.5 * node_values_derivative_[1] +
                1.5 * (node_values_[1] - node_values_[0]) /
                (nodes_[1] - nodes_[0]);
            break;

        case SplineBoundaryCondition::naturalZeroDerivative:
            throw AmiException("Natural boundary condition with zero
                "derivative is not allowed for Hermite splines.");

        case SplineBoundaryCondition::periodic:
            if (node_derivative_by_FD_)
                node_values_derivative_[0] =
                    (node_values_[1] - node_values_[last - 1]) /
                    (nodes_[1] - nodes_[0] + nodes_[last] - nodes_[last - 1]);
            break;
    }

    /* ...and the last node (1-sided FD). */
    switch (lastNodeBC_) {
        case SplineBoundaryCondition::given:
            if (node_derivative_by_FD_)
                /* 1-sided FD */
                node_values_derivative_[last] =
                    (node_values_[last] - node_values_[last - 1]) /
                    (nodes_[last] - nodes_[last - 1]);
            break;

        case SplineBoundaryCondition::zeroDerivative:
            node_values_derivative_[0] = 0;
            break;

        case SplineBoundaryCondition::natural:
            node_values_derivative_[last] =
                -0.5 * node_values_derivative_[last - 1] +
                1.5 * (node_values_[last] - node_values_[last - 1]) /
                (nodes_[last] - nodes_[last - 1]);
            break;

        case SplineBoundaryCondition::naturalZeroDerivative:
            throw AmiException("Natural boundary condition with zero
                "derivative is not allowed for Hermite splines.");

        case SplineBoundaryCondition::periodic:
            if (node_derivative_by_FD_)
                node_values_derivative_[last] =
                    (node_values_[1] - node_values_[last - 1]) /
                    (nodes_[1] - nodes_[0] + nodes_[last] - nodes_[last - 1]);
            break;
    }
}

void HermiteSpline::computeCoefficients() {
    /* Allocate space for the coefficients for Horner's method.
     * They are stored in the vector as
     * [d_0, c_0, b_0, a_0, d_1, c_1, ... , b_{n_nodes-1}, a_{n_nodes-1}] */
    coefficients.resize(4 * (n_nodes() - 1), 0.0);
    realtype len;

    /* Compute the coefficients of the spline polynomials:
     * spline(t) = a * t**3 + b * t**2 + c * t + d
     *           = d + t * (c + t * (b + t * a))
     * with coefficients[4 * i_node + (0, 1, 2, 3)] = (d, c, b, a)
     * */

    for (int i_node = 0; i_node < n_nodes() - 1; i_node++) {
        /* Get the length of the interval. Yes, we could save computation time
         * by exploiting equidistant spacing, but we're talking about <1k FLOPs
         * for sure, no matter what model. Screw it. */
        len = nodes_[i_node + 1] - nodes_[i_node];

        /* Coefficients for cubic Hermite polynomials */
        coefficients[4 * i_node] = node_values_[i_node];
        coefficients[4 * i_node + 1] = len*node_values_derivative_[i_node];
        coefficients[4 * i_node + 2] = - 3 * node_values_[i_node]
            - 2 * len * node_values_derivative_[i_node]
            + 3 * node_values_[i_node + 1]
            - len * node_values_derivative_[i_node + 1];
        coefficients[4 * i_node + 3] = 2 * node_values_[i_node]
            + len * node_values_derivative_[i_node]
            - 2 * node_values_[i_node + 1]
            + len * node_values_derivative_[i_node + 1];
    }

    /* Take care of coefficients for extrapolation */
    computeCoefficientsExtrapolation();
}

void HermiteSpline::computeCoefficientsExtrapolation() {
    /* Do we want to extrapolate at all? */
    bool needExtrapolationCoefficients =
        (firstNodeEP_ == SplineExtrapolation::constant ||
         firstNodeEP_ == SplineExtrapolation::linear) &&
        (lastNodeEP_ == SplineExtrapolation::constant ||
         lastNodeEP_ == SplineExtrapolation::linear);
    if (!needExtrapolationCoefficients)
        return;

    coefficients_extrapolate.resize(4, 0.0);

    int last = n_nodes() - 1;

    /* Beyond the spline nodes, we need to extrapolate using a * t + b.
     * Those coefficients are stored as [b_first, a_first, b_last, a_last] */
    switch (firstNodeEP_) {
        case SplineExtrapolation::constant:
            coefficients_extrapolate[0] = node_values_[0];
            coefficients_extrapolate[1] = 0;
            break;

        case SplineExtrapolation::linear:
            coefficients_extrapolate[0] = node_values_[0]
                - nodes_[0] * node_values_derivative_[0];
            coefficients_extrapolate[1] = node_values_derivative_[0];
            break;

        default:
            /* We don't need specific coefficients in the cases of:
             * noExtrapolation, polynomial, periodic*/
            break;
    }
    switch (lastNodeEP_) {
        case SplineExtrapolation::constant:
            coefficients_extrapolate[2] = node_values_[last];
            coefficients_extrapolate[3] = 0;
            break;

        case SplineExtrapolation::linear:
            coefficients_extrapolate[2] = node_values_[last]
                - nodes_[last] * node_values_derivative_[last];
            coefficients_extrapolate[3] = node_values_derivative_[last];
            break;

        default:
            /* We don't need specific coefficients in the cases of:
             * noExtrapolation, polynomial, periodic*/
            break;
    }
)

void HermiteSpline::computeCoefficientsSensi(int nplist, int spline_offset,
                                             realtype *dspline_valuesdp,
                                             realtype *dspline_slopesdp) {
    /* Allocate space for the coefficients *
     * They are stored in the vector as
     * [ D[d_0, p0], D[c_0, p0], D[b_0, p0], D[a_0, p0], D[d_1, p0],
     *   ... ,
     *   D[b_{n_nodes-1}, p0], D[a_{n_nodes-1}, p0],
     *   D[d_0, p1], D[c_0, p1], ...
     *   ..., D[b_{n_nodes-1}, p{nplist-1}, D[a_{n_nodes-1}, p{nplist-1}]
     * ] */
    int n_spline_coefficients = 4 * (n_nodes() - 1);
    coefficients_sensi.resize(n_spline_coefficients * nplist, 0.0);

    /**
      * We're using short hand notation for some node values or slopes, based on
      * the notation used on https://en.wikipedia.org/wiki/Cubic_Hermite_spline
      * In brief: "p" denotes the current (k-th) spline node value,
      * "m" its tangent or slope, "s" in front the sensitivity, "1" at the end
      * means the following node (" + 1"), so "smk1" is the sensitivity of the
      * slope at node k + 1, w.r.t. to the current parameter (looping index).
      * */
    realtype len, len_m, len_p;

    /* Parametric derivatives of splines are splines again.
     * We compute the coefficients for those polynomials now. */
    for (int i_node = 0; i_node < n_nodes() - 2; i_node++) {
        /* Get the length of the interval. */
        len = nodes_[i_node + 1] - nodes_[i_node];
        len_m = (i_node > 0) ? nodes_[i_node + 1] - nodes_[i_node - 1]
                             : len;
        len_p = (i_node < n_nodes() - 2) ? nodes_[i_node + 2] - nodes_[i_node]
                                         : len;

        /* As computing the coefficient is a mess, it's in another function */
        for (int ip = 0; ip < nplist; ip++)
            getCoeffsSensiLowlevel(ip, i_node, nplist, n_spline_coefficients,
                                   spline_offset, len, len_m, len_p,
                                   dspline_valuesdp, dspline_slopesdp,
                                   coefficients_sensi.data(),
                                   coefficients_extrapolate_sensi.data());
    }

    /* We need the coefficients for extrapolating beyond the spline domain */
    computeCoefficientsExtrapolationSensi(nplist, spline_offset,
                                          dspline_valuesdp, dspline_slopesdp);
}

void HermiteSpline::computeCoefficientsExtrapolationSensi(int nplist,
    int spline_offset, realtype *dspline_valuesdp, realtype *dspline_slopesdp) {

    /* Do we want to extrapolate at all? */
    bool needExtrapolationCoefficients =
        (firstNodeEP_ == SplineExtrapolation::constant ||
         firstNodeEP_ == SplineExtrapolation::linear) &&
        (lastNodeEP_ == SplineExtrapolation::constant ||
         lastNodeEP_ == SplineExtrapolation::linear);
    if (!needExtrapolationCoefficients)
        return;

    /* Beyond the spline nodes, we need to extrapolate using a * t + b.
     * Those coefficients are stored as
     * [ D[b_first, p0], D[a_first, p0], D[b_last, p0], D[a_last, p0],
     *   D[b_first, p1], ... D[a_last, p{nplist-1}]
     * ] */
    coefficients_extrapolate_sensi.resize(4 * nplist, 0.0);

    realtype sp0;
    realtype sm0;
    for (int ip = 0; ip < nplist; ip++) {
        sp0 = dspline_valuesdp[spline_offset + ip];
        switch (firstNodeEP_) {
            /* This whole switch-case-if-else-if-thing could be moved
             * outside the loop, I know. Yet, it's at most some thousand
             * if's done once in the program for saving many lines of code
             * and getting a much clearer code structure. */
            case SplineExtrapolation::constant:
                sm0 = 0;
                break;

            case SplineExtrapolation::linear:
                if (get_node_derivative_by_FD() &&
                    firstNodeBC_ == SplineBoundaryCondition::given) {
                    sm0 = (dspline_valuesdp[spline_offset + ip + nplist] - sp0)
                        / (nodes_[1] - nodes_[0]);

                } else if (get_node_derivative_by_FD() &&
                           firstNodeBC_ == SplineBoundaryCondition::natural) {
                    throw AmiException("Natural boundary condition for "
                        "Hermite splines with linear extrapolation is "
                        "not yet implemented.");

                } else if (!get_node_derivative_by_FD() &&
                           firstNodeBC_ == SplineBoundaryCondition::given) {
                    sm0 = dspline_slopesdp[spline_offset + ip];
                    throw AmiException("Natural boundary condition for "
                        "Hermite splines with linear extrapolation is "
                        "not yet implemented.");

                } else if (!get_node_derivative_by_FD() &&
                           firstNodeBC_ == SplineBoundaryCondition::natural) {

                } else {
                    throw AmiException("Some weird combination of spline boundary "
                        "condition, extrapolation and finite differecnces was "
                        "passed which should bo ne allowed.");
                }
                break;

            default:
                /* We don't need specific coefficients in the cases of:
                 * noExtrapolation, polynomial, periodic*/
                break;
        }
        /* Write them to the vector */
        coefficients_extrapolate_sensi[4 * ip] = sp0 - sm0 * nodes_[0];
        coefficients_extrapolate_sensi[4 * ip + 1] = sm0;
    }

    realtype sp_end;
    realtype sm_end;
    for (int ip = 0; ip < nplist; ip++) {
        sp_end = dspline_valuesdp[spline_offset + ip + (n_nodes()-1) * nplist];
        switch (lastNodeEP_) {
            /* This whole switch-case-if-else-if-thing could be moved
             * outside the loop, I know. Yet, it's at most some thousand
             * if's done once in the program for saving many lines of code
             * and getting a much clearer code structure. */
            case SplineExtrapolation::constant:
                sm_end = 0;
                break;

            case SplineExtrapolation::linear:
                if (get_node_derivative_by_FD() &&
                    lastNodeBC_ == SplineBoundaryCondition::given) {
                    sm_end = (sp_end - dspline_valuesdp[spline_offset +
                              ip + (n_nodes() - 2) * nplist])
                             / (nodes_[n_nodes() - 1] - nodes_[n_nodes() - 2]);

                } else if (get_node_derivative_by_FD() &&
                           lastNodeBC_ == SplineBoundaryCondition::natural) {
                    throw AmiException("Natural boundary condition for "
                        "Hermite splines with linear extrapolation is "
                        "not yet implemented.");

                } else if (!get_node_derivative_by_FD() &&
                           lastNodeBC_ == SplineBoundaryCondition::given) {
                    sm_end = dspline_slopesdp[spline_offset + ip +
                             (n_nodes() - 1) * nplist];

                } else if (!get_node_derivative_by_FD() &&
                           lastNodeBC_ == SplineBoundaryCondition::natural) {
                    throw AmiException("Natural boundary condition for "
                        "Hermite splines with linear extrapolation is "
                        "not yet implemented.");

                } else {
                    throw AmiException("Some weird combination of spline boundary "
                        "condition, extrapolation and finite differecnces was "
                        "passed which should bo ne allowed.");
                }
                break;

            default:
                /* We don't need specific coefficients in the cases of:
                 * noExtrapolation, polynomial, periodic*/
                break;
        }
        /* Write them to the vector */
        coefficients_extrapolate_sensi[4 * ip + 2] = sp_end -
            sm_end * nodes_[n_nodes() - 1];
        coefficients_extrapolate_sensi[4 * ip + 3] = sm_end;
    }
}

void HermiteSpline::getCoeffsSensiLowlevel(int ip, int i_node, int nplist, int n_spline_coefficients,
                                           int spline_offset, realtype len, realtype len_m,
                                           realtype len_p, realtype *dnodesdp,
                                           realtype *dslopesdp, realtype *coeffs,
                                           realtype *coeffs_extrapol) {
    /**
     * We're using the short hand notation for node values and slopes from
     * computeCoefficientsSensi() here. See this function for documentation.
     * */
    int node_offset = spline_offset + ip;
    int last = n_nodes() - 1;
    double spk = dnodesdp[node_offset + i_node * nplist];
    double spk1 = dnodesdp[node_offset + (i_node + 1) * nplist];
    double smk;
    double smk1;

    /* Get sensitivities of slopes. Depending on finite differences are used
     * or not, this may be a bit cumbersome now... */
    if (get_node_derivative_by_FD()) {
        if (i_node == 0) {
            /* Are we at the fist node? What's the boundary condition? */
            if (firstNodeBC_ == SplineBoundaryCondition::zeroDerivative) {
                smk = 0;
            } else if (firstNodeBC_ == SplineBoundaryCondition::given) {
                smk = (spk1 - spk) / len_m;
            } else if (firstNodeBC_ == SplineBoundaryCondition::natural) {
                throw AmiException("Natural boundary condition for Hermite "
                   "splines is not yet implemented.");
            } else if (firstNodeBC_ == SplineBoundaryCondition::periodic) {
                smk = (spk1 - dnodesdp[node_offset + last * nplist]) /
                      (len_m + nodes_[last] - nodes_[last - 1]);
            } else {
                /* must be SplineBoundaryCondition::naturalZeroDerivative*/
                throw AmiException("Natural boundary condition with zero "
                    "derivative is prohibited for Hermite splines.");
            }
            smk1 = (dnodesdp[node_offset + (i_node + 2) * nplist] - spk) / len_p;

        } else if (i_node == n_nodes() - 2) {
            /* Are we at the last node? What's the boundary condition? */
            smk = (spk1 - dnodesdp[node_offset + (i_node - 1) * nplist]) / len_m;
            if (lastNodeBC_ == SplineBoundaryCondition::zeroDerivative) {
                smk1 = 0;
            } else if (lastNodeBC_ == SplineBoundaryCondition::given) {
                smk1 = (spk1 - spk) / len_p;
            } else if (lastNodeBC_ == SplineBoundaryCondition::natural) {
                throw AmiException("Natural boundary condition for Hermite "
                   "splines is not yet implemented.");
            } else if (lastNodeBC_ == SplineBoundaryCondition::periodic) {
                smk1 = (spk1 - dnodesdp[node_offset]) /
                       (len_p + nodes_[1] - nodes_[0]);
            } else {
                /* must be SplineBoundaryCondition::naturalZeroDerivative*/
                throw AmiException("Natural boundary condition with zero "
                    "derivative is prohibited for Hermite splines.");
            }
            smk1 = (dnodesdp[node_offset + (i_node + 2) * nplist] - spk) / len_p;

        } else {
            /* We're somewhere in between. That's fine. */
            smk = (spk1 - dnodesdp[node_offset + (i_node - 1) * nplist]) / len_m;
            smk1 = (dnodesdp[node_offset + (i_node + 2) * nplist] - spk) / len_p;
        }
    } else {
        /* The slopes are explicitly given, easiest case... */
        smk = dslopesdp[node_offset + i_node * nplist];
        smk1 = dslopesdp[node_offset + (i_node + 1) * nplist];

        /* For the nodes at the boundary, we have to take care of the bc */
        if (i_node == 0 &&
            firstNodeBC_ == SplineBoundaryCondition::zeroDerivative)
            smk = 0;
        if (i_node == n_nodes() - 2 &&
            lastNodeBC_ == SplineBoundaryCondition::zeroDerivative)
            smk1 = 0;
    }

    /* Compute the actual coefficients */
    coeffs[ip * n_spline_coefficients + 4 * i_node] = spk;
    coeffs[ip * n_spline_coefficients + 4 * i_node + 1] = len * smk;
    coeffs[ip * n_spline_coefficients + 4 * i_node + 2] = 3 * (spk1  - spk) - len * (2 * smk + smk1);
    coeffs[ip * n_spline_coefficients + 4 * i_node + 3] = 2 * (spk - spk1) + len * (smk + smk1);
}

void HermiteSpline::computeFinalValue() {
    /* We need to compute the final value of the spline, depending on its
     * boundary condition and the extrapolation option. */
    realtype finalValue;
    if ((lastNodeEP_ == SplineExtrapolation::constant) ||
        (lastNodeBC_ == SplineBoundaryCondition::zeroDerivative &&
         lastNodeEP_ == SplineExtrapolation::linear)) {
        finalValue = coefficients_extrapolate[2];
    } else if (lastNodeEP_ == SplineExtrapolation::linear) {
        if (coefficients_extrapolate[2] < 0) {
            finalValue = -INFINITY;
        } else if (coefficients_extrapolate[2] > 0) {
            finalValue = INFINITY;
        } else {
            finalValue = coefficients_extrapolate[2];
        }
    } else if (lastNodeEP_ == SplineExtrapolation::polynomial) {
        int last = 4 * (n_nodes() - 1) - 1;
        if (coefficients[last] < 0) {
            finalValue = -INFINITY;
        } else if (coefficients[last] > 0) {
            finalValue = INFINITY;
        } else {
            if (coefficients[last - 1] < 0) {
                finalValue = -INFINITY;
            } else if (coefficients[last - 1] > 0) {
                finalValue = INFINITY;
            } else {
                if (coefficients[last - 2] < 0) {
                    finalValue = -INFINITY;
                } else if (coefficients[last - 2] > 0) {
                    finalValue = INFINITY;
                } else {
                    finalValue = coefficients[last - 3];
                }
            }
        }
    } else {
        /* Periodic: will not yield a steady state */
        finalValue = NAN;
    }
    setFinalValue(finalValue)
}

void HermiteSpline::computeFinalSensitivity(int nplist, int spline_offset,
                                            realtype *dspline_valuesdp,
                                            realtype *dspline_slopesdp) {
    /* We need to compute the final value of the spline, depending on its
     * boundary condition and the extrapolation option. */
    std::vector<realtype> finalSensitivity(nplist, 0);
    if ((lastNodeEP_ == SplineExtrapolation::constant) ||
        (lastNodeBC_ == SplineBoundaryCondition::zeroDerivative &&
         lastNodeEP_ == SplineExtrapolation::linear)) {
        for (int ip = 0; ip < nplist; ip++)
            // TODO:
        finalValue = coefficients_extrapolate[2];
    } else if (lastNodeEP_ == SplineExtrapolation::linear) {
        if (coefficients_extrapolate[2] < 0) {
            finalValue = -INFINITY;
        } else if (coefficients_extrapolate[2] > 0) {
            finalValue = INFINITY;
        } else {
            finalValue = coefficients_extrapolate[2];
        }
    } else if (lastNodeEP_ == SplineExtrapolation::polynomial) {
        /* Yes, that's not correct. But I don't see any good reason for
         * implementing a case, which anybody with more than a dead fish
         * between the ears will never use. */
         std::fill(finalSensitivity.begin(), finalSensitivity.end(), NAN);
    } else {
        /* Periodic: will not yield a steady state */
        std::fill(finalSensitivity.begin(), finalSensitivity.end(), NAN);
    }
}

realtype HermiteSpline::getValue(const double t) {
    /* Is this a steady state computation? */
    if std::isinf(t)
        return finalValue;

    /* Compute the spline value */
    int i_node = 0;
    realtype len = nodes_[1] - nodes_[0];

    /* Are we past the last node? Extrapolate! */
    if (t > nodes_[n_nodes() - 1]) {
        switch (lastNodeEP_) {
            case SplineExtrapolation:noExtrapolation:
                throw AmiException("Trying to evaluate spline after last "
                    "spline node, but spline has been specified to not allow "
                    "extrapolation.");

            case SplineExtrapolation::constant:
                return coefficients_extrapolate[2];

            case SplineExtrapolation::linear:
                return coefficients_extrapolate[2] +
                    t * coefficients_extrapolate[3];

            case SplineExtrapolation::polynomial:
                /* Evaluate last interpolation polynomial */
                i_node = n_nodes() - 2;
                len = nodes_[i_node + 1] - nodes_[i_node];
                return evaluatePolynomial((t - nodes_[i_node]) / len,
                                          &(coefficients[i_node * 4]));

            case SplineExtrapolation::periodic:
                len = nodes_[n_nodes() - 1] - nodes_[0];
                return getValue(t - len * std::floor((t - nodes_[0]) /
                    len));
        }
    }

    /* Are we before the first node? Extrapolate! */
    if (t < nodes_[0]) {
        switch (firstNodeEP_) {
            case SplineExtrapolation:noExtrapolation:
                throw AmiException("Trying to evaluate spline before first "
                    "spline node, but spline has been specified to not allow "
                    "extrapolation.");

            case SplineExtrapolation::constant:
                return coefficients_extrapolate[0];

            case SplineExtrapolation::linear:
                return coefficients_extrapolate[0] +
                    t * coefficients_extrapolate[1];

            case SplineExtrapolation::polynomial:
                /* Evaluate last interpolation polynomial */
                len = nodes_[1] - nodes_[0];
                return evaluatePolynomial((t - nodes_[0]) / len,
                                          &(coefficients[0]));

            case SplineExtrapolation::periodic:
                return getValue(t + len + len * std::floor((t - nodes_[0]) /
                    len));
        }

    }

    /* Get the spline interval which we need */
    if (get_equidistant_spacing()) {
        /* equidistant spacing: just compute the interval */
        i_node = std::floor((t - nodes_[0]) / len);
    } else {
        /* no equidistant spacing: we need to iterate */
        while (nodes_[i_node + 1] < t) {
            i_node++;
        }
        len = nodes_[i_node + 1] - nodes_[i_node];
    }

    /* Evaluate the interpolation polynomial */
    return evaluatePolynomial((t - nodes_[i_node]) / len,
                              &(coefficients[i_node * 4]));
}

realtype HermiteSpline::getSensitivity(const double t, const int ip) {
    /* Is this a steady state computation? */
    if std::isinf(t)
        return finalSensitivity[ip];

    /* Compute the spline value */
    int i_node = 0;
    realtype len;

    /* Compute the parametric derivative of the spline value */
    if (t > nodes_[n_nodes() - 1]) {
        /* Are we past the last node? Extrapolate! */
        switch (lastNodeEP_) {
            case SplineExtrapolation:noExtrapolation:
                throw AmiException("Trying to evaluate spline sensitivity "
                    "after last spline node, but spline has been specified "
                    "to not allow extrapolation.");

            case SplineExtrapolation::constant:
                return coefficients_extrapolate_sensi[4 * ip + 2];

            case SplineExtrapolation::linear:
                return coefficients_extrapolate_sensi[4 * ip + 2]
                    + t * coefficients_extrapolate_sensi[4 * ip + 3];

            case SplineExtrapolation::polynomial:
                /* Evaluate last interpolation polynomial */
                i_node = n_nodes() - 2;
                len = nodes_[i_node + 1] - nodes_[i_node];
                return evaluatePolynomial((t - nodes_[i_node]) / len,
                                  &(coefficients_sensi[ip * (n_nodes() - 1) * 4
                                                       + i_node * 4]));

            case SplineExtrapolation::periodic:
                len = nodes_[n_nodes() - 1] - nodes_[0];
                return getSensitivity(t - len * std::floor((t - nodes_[0]) /
                    len), ip);
        }

    } else if (t < nodes_[0]) {
        /* Are we before the first node? Extrapolate! */
        switch (firstNodeEP_) {
            case SplineExtrapolation:noExtrapolation:
                throw AmiException("Trying to evaluate spline before first "
                    "spline node, but spline has been specified to not allow "
                    "extrapolation.");

            case SplineExtrapolation::constant:
                return coefficients_extrapolate_sensi[4 * ip + 0];

            case SplineExtrapolation::linear:
                return coefficients_extrapolate_sensi[4 * ip + 0]
                    + t * coefficients_extrapolate_sensi[4 * ip + 1];

            case SplineExtrapolation::polynomial:
                /* Evaluate last interpolation polynomial */
                len = nodes_[1] - nodes_[0];
                return evaluatePolynomial((t - nodes_[0]) / len,
                    &(coefficients_sensi[ip * (n_nodes() - 1) * 4]));

            case SplineExtrapolation::periodic:
                return getSensitivity(t + len * (1 + std::floor((t - nodes_[0])
                    / len)), ip);
        }

    } else {
        /* Get the spline interval which we need */
        if (get_equidistant_spacing()) {
            /* equidistant spacing: just compute the interval */
            len = nodes_[1] - nodes_[0];
            i_node = std::floor((t - nodes_[0]) / len);
        } else {
            /* no equidistant spacing: we need to iterate */
            while (nodes_[i_node + 1] < t)
                i_node++;

            len = nodes_[i_node + 1] - nodes_[i_node];
        }

        return evaluatePolynomial((t - nodes_[i_node]) / len,
            &(coefficients_sensi[ip * (n_nodes() - 1) * 4 + i_node * 4]));
    }

}

} // namespace amici
