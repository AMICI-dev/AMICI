#include "amici/splinefunctions.h"
#include "amici/amici.h"
#include "amici/defines.h"
#include "amici/exception.h"
#include "amici/vector.h"

#include <algorithm> // std::min
#include <cmath>
#include <vector>

namespace amici {

static realtype
evaluate_polynomial(realtype const x, gsl::span<const realtype> coeff)
{
    /* Use Horner's method (https://en.wikipedia.org/wiki/Horner%27s_method)
     * for numerical efficiency:
     *
     * spline(t) = a * t**3 + b * t**2 + c * t + d
     *           = d + t * (c + t * (b + t * a))
     * with coeff[0, 1, 2, 3] = [d, c, b, a]
     */
    assert(coeff.size() >= 4);
    auto coeff_p = coeff.data();
    return coeff_p[0] + x * (coeff_p[1] + x * (coeff_p[2] + x * coeff_p[3]));
}

AbstractSpline::AbstractSpline(std::vector<realtype> nodes,
                               std::vector<realtype> node_values,
                               bool equidistant_spacing,
                               bool logarithmic_parametrization)
  : nodes_(std::move(nodes))
  , node_values_(std::move(node_values))
  , equidistant_spacing_(equidistant_spacing)
  , logarithmic_parametrization_(logarithmic_parametrization)
{

    /* we want to set the number of nodes */
    auto n_nodes_ = static_cast<int>(node_values_.size());

    /* In case we have equidistant spacing, compute node locations */
    if (equidistant_spacing_) {
        if (nodes_.size() != 2)
            throw AmiException("Splines with equidistant spacing need a nodes "
                               "vector with two elements (first/last node).");
        realtype node_start = nodes_[0];
        realtype node_step = (nodes_[1] - nodes_[0]) / (n_nodes_ - 1);
        nodes_.resize(n_nodes_);
        nodes_[n_nodes_ - 1] = nodes_[1];
        for (int i_node = 0; i_node < n_nodes_ - 1; i_node++)
            nodes_[i_node] = node_start + i_node * node_step;
    }
    else if (nodes_.size() != node_values_.size()) {
        throw std::invalid_argument(
            "Number of nodes and number of node_values do not match.");
    }

    if (logarithmic_parametrization_) {
        for (auto& node_value : node_values_)
            node_value = std::log(node_value);
    }
}

realtype
AbstractSpline::get_final_value() const
{
    return final_value_;
}

void
AbstractSpline::set_final_value(realtype finalValue)
{
    final_value_ = finalValue;
}

realtype
AbstractSpline::get_final_sensitivity(const int ip) const
{
    return final_sensitivity_[ip];
}

void
AbstractSpline::set_final_sensitivity(std::vector<realtype> finalSensitivity)
{
    final_sensitivity_ = std::move(finalSensitivity);
}

bool
AbstractSpline::get_equidistant_spacing() const
{
    return equidistant_spacing_;
}

void
AbstractSpline::set_equidistant_spacing(bool equidistant_spacing)
{
    equidistant_spacing_ = equidistant_spacing;
}

bool
AbstractSpline::get_logarithmic_parametrization() const
{
    return logarithmic_parametrization_;
}

void
AbstractSpline::set_logarithmic_parametrization(
  bool logarithmic_parametrization)
{
    logarithmic_parametrization_ = logarithmic_parametrization;
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
                             bool logarithmic_parametrization)
  : AbstractSpline(std::move(nodes),
                   std::move(node_values),
                   equidistant_spacing,
                   logarithmic_parametrization)
  , node_values_derivative_(std::move(node_values_derivative))
  , first_node_bc_(firstNodeBC)
  , last_node_bc_(lastNodeBC)
  , first_node_ep_(firstNodeExtrapol)
  , last_node_ep_(lastNodeExtrapol)
  , node_derivative_by_FD_(node_derivative_by_FD)
{
    if (not node_derivative_by_FD_ and
        node_values_derivative_.size() != nodes_.size()) {
        throw std::invalid_argument(
          "Size of node_values_derivative does not match number of nodes.");
    }

    /* We may have to compute the derivatives at the nodes */
    handle_inner_derivatives();
    /* First and last node need to be handled separately */
    handle_boundary_conditions();
}

void
HermiteSpline::handle_inner_derivatives()
{
    /* If values of the derivative at the nodes are to be computed by finite
     * differences, we have to fill up node_values_derivative_ */
    if (node_derivative_by_FD_) {
        node_values_derivative_.resize(n_nodes(), 0.0);
        if (equidistant_spacing_) {
            realtype hx2 = 2 * (nodes_[1] - nodes_[0]);
            for (int i_node = 1; i_node < n_nodes() - 1; i_node++)
              node_values_derivative_[i_node] =
                (node_values_[i_node + 1] - node_values_[i_node - 1]) / hx2;
        } else {
            for (int i_node = 1; i_node < n_nodes() - 1; i_node++) {
              realtype dleft =
                (node_values_[i_node] - node_values_[i_node - 1]) /
                (nodes_[i_node] - nodes_[i_node - 1]);
              realtype dright =
                (node_values_[i_node + 1] - node_values_[i_node]) /
                (nodes_[i_node + 1] - nodes_[i_node]);
              node_values_derivative_[i_node] = (dleft + dright) / 2;
            }
        }
    }
}

void
HermiteSpline::handle_boundary_conditions()
{
    int last = n_nodes() - 1;

    /* We have to take special care of the first node */
    switch (first_node_bc_) {
        case SplineBoundaryCondition::given:
            if (node_derivative_by_FD_)
                /* 1-sided FD */
                node_values_derivative_[0] =
                  (node_values_[1] - node_values_[0]) / (nodes_[1] - nodes_[0]);
            break;

        case SplineBoundaryCondition::zeroDerivative:
            node_values_derivative_[0] = 0;
            break;

        case SplineBoundaryCondition::natural:
            node_values_derivative_[0] = -0.5 * node_values_derivative_[1] +
                                         1.5 *
                                           (node_values_[1] - node_values_[0]) /
                                           (nodes_[1] - nodes_[0]);
            break;

        case SplineBoundaryCondition::naturalZeroDerivative:
            throw AmiException(
              "Natural boundary condition with zero "
              "derivative is not allowed for Hermite splines.");

        case SplineBoundaryCondition::periodic:
            if (node_derivative_by_FD_)
                node_values_derivative_[0] =
                  (node_values_[1] - node_values_[last - 1]) /
                  (nodes_[1] - nodes_[0] + nodes_[last] - nodes_[last - 1]);
            break;
        default:
            throw AmiException("Invalid value for boundary condition.");
    }

    /* ...and the last node (1-sided FD). */
    switch (last_node_bc_) {
        case SplineBoundaryCondition::given:
            if (node_derivative_by_FD_)
                /* 1-sided FD */
                node_values_derivative_[last] =
                  (node_values_[last] - node_values_[last - 1]) /
                  (nodes_[last] - nodes_[last - 1]);
            break;

        case SplineBoundaryCondition::zeroDerivative:
            node_values_derivative_[last] = 0;
            break;

        case SplineBoundaryCondition::natural:
            node_values_derivative_[last] =
              -0.5 * node_values_derivative_[last - 1] +
              1.5 * (node_values_[last] - node_values_[last - 1]) /
                (nodes_[last] - nodes_[last - 1]);
            break;

        case SplineBoundaryCondition::naturalZeroDerivative:
            throw AmiException(
              "Natural boundary condition with zero "
              "derivative is not allowed for Hermite splines.");

        case SplineBoundaryCondition::periodic:
            if (node_derivative_by_FD_)
                node_values_derivative_[last] =
                  (node_values_[1] - node_values_[last - 1]) /
                  (nodes_[1] - nodes_[0] + nodes_[last] - nodes_[last - 1]);
            break;
        default:
            throw AmiException("Invalid value for boundary condition.");
    }
}

void
HermiteSpline::compute_coefficients()
{
    /* Allocate space for the coefficients for Horner's method.
     * They are stored in the vector as
     * [d_0, c_0, b_0, a_0, d_1, c_1, ... , b_{n_nodes-1}, a_{n_nodes-1}] */
    coefficients.resize(4 * (n_nodes() - 1), 0.0);

    /* Compute the coefficients of the spline polynomials:
     * spline(t) = a * t**3 + b * t**2 + c * t + d
     *           = d + t * (c + t * (b + t * a))
     * with coefficients[4 * i_node + (0, 1, 2, 3)] = (d, c, b, a)
     */

    for (int i_node = 0; i_node < n_nodes() - 1; i_node++) {
        /* Get the length of the interval. Yes, we could save computation time
         * by exploiting equidistant spacing, but we're talking about <1k FLOPs
         * for sure, no matter what model. Screw it. */
        realtype len = nodes_[i_node + 1] - nodes_[i_node];

        /* Coefficients for cubic Hermite polynomials */
        coefficients[4 * i_node] = node_values_[i_node];
        coefficients[4 * i_node + 1] = len * node_values_derivative_[i_node];
        coefficients[4 * i_node + 2] =
          -3 * node_values_[i_node] -
          2 * len * node_values_derivative_[i_node] +
          3 * node_values_[i_node + 1] -
          len * node_values_derivative_[i_node + 1];
        coefficients[4 * i_node + 3] =
          2 * node_values_[i_node] + len * node_values_derivative_[i_node] -
          2 * node_values_[i_node + 1] +
          len * node_values_derivative_[i_node + 1];
    }

    /* Take care of coefficients for extrapolation */
    compute_coefficients_extrapolation();
}

void
HermiteSpline::compute_coefficients_extrapolation()
{
    /* Do we want to extrapolate at all? */
    bool needExtrapolationCoefficients =
      first_node_ep_ == SplineExtrapolation::constant ||
      first_node_ep_ == SplineExtrapolation::linear ||
      last_node_ep_ == SplineExtrapolation::constant ||
      last_node_ep_ == SplineExtrapolation::linear;
    if (!needExtrapolationCoefficients)
        return;

    coefficients_extrapolate.resize(4, 0.0);

    int last = n_nodes() - 1;

    /* Beyond the spline nodes, we need to extrapolate using a * t + b.
     * Those coefficients are stored as [b_first, a_first, b_last, a_last] */
    switch (first_node_ep_) {
        case SplineExtrapolation::constant:
            coefficients_extrapolate[0] = node_values_[0];
            coefficients_extrapolate[1] = 0;
            break;

        case SplineExtrapolation::linear:
            coefficients_extrapolate[0] =
              node_values_[0] - nodes_[0] * node_values_derivative_[0];
            coefficients_extrapolate[1] = node_values_derivative_[0];
            break;

        default:
            /* We don't need specific coefficients in the cases of:
             * noExtrapolation, polynomial, periodic*/
            break;
    }
    switch (last_node_ep_) {
        case SplineExtrapolation::constant:
            coefficients_extrapolate[2] = node_values_[last];
            coefficients_extrapolate[3] = 0;
            break;

        case SplineExtrapolation::linear:
            coefficients_extrapolate[2] =
              node_values_[last] - nodes_[last] * node_values_derivative_[last];
            coefficients_extrapolate[3] = node_values_derivative_[last];
            break;

        default:
            /* We don't need specific coefficients in the cases of:
             * noExtrapolation, polynomial, periodic*/
            break;
    }
}

void
HermiteSpline::compute_coefficients_sensi(int nplist,
                                          int spline_offset,
                                          gsl::span<realtype> dspline_valuesdp,
                                          gsl::span<realtype> dspline_slopesdp)
{
    /*
     * Allocate space for the coefficients
     * They are stored in the vector as
     * [ D[d_0, p0], D[c_0, p0], D[b_0, p0], D[a_0, p0], D[d_1, p0],
     *   ... ,
     *   D[b_{n_nodes-1}, p0], D[a_{n_nodes-1}, p0],
     *   D[d_0, p1], D[c_0, p1], ...
     *   ..., D[b_{n_nodes-1}, p{nplist-1}, D[a_{n_nodes-1}, p{nplist-1}]
     * ]
     */
    int n_spline_coefficients = 4 * (n_nodes() - 1);
    coefficients_sensi.resize(n_spline_coefficients * nplist, 0.0);

    /*
     * We're using short hand notation for some node values or slopes, based on
     * the notation used on https://en.wikipedia.org/wiki/Cubic_Hermite_spline
     * In brief: "p" denotes the current (k-th) spline node value,
     * "m" its tangent or slope, "s" in front the sensitivity, "1" at the end
     * means the following node (" + 1"), so "smk1" is the sensitivity of the
     * slope at node k + 1, w.r.t. to the current parameter (looping index).
     */

    /* Parametric derivatives of splines are splines again.
     * We compute the coefficients for those polynomials now. */
    for (int i_node = 0; i_node < n_nodes() - 1; i_node++) {
        /* Get the length of the interval. */
        realtype len = nodes_[i_node + 1] - nodes_[i_node];
        realtype len_m =
          (i_node > 0) ? nodes_[i_node + 1] - nodes_[i_node - 1] : len;
        realtype len_p =
          (i_node < n_nodes() - 2) ? nodes_[i_node + 2] - nodes_[i_node] : len;

        /* As computing the coefficient is a mess, it's in another function */
        for (int ip = 0; ip < nplist; ip++)
            get_coeffs_sensi_lowlevel(ip,
                                      i_node,
                                      nplist,
                                      n_spline_coefficients,
                                      spline_offset,
                                      len,
                                      len_m,
                                      len_p,
                                      dspline_valuesdp,
                                      dspline_slopesdp,
                                      coefficients_sensi);
    }

    /* We need the coefficients for extrapolating beyond the spline domain */
    compute_coefficients_extrapolation_sensi(
      nplist, spline_offset, dspline_valuesdp, dspline_slopesdp);
}

void
HermiteSpline::compute_coefficients_extrapolation_sensi(
  int nplist,
  int spline_offset,
  gsl::span<realtype> dspline_valuesdp,
  gsl::span<realtype> dspline_slopesdp)
{

    /* Do we want to extrapolate at all? */
    bool needExtrapolationCoefficients =
      first_node_ep_ == SplineExtrapolation::constant ||
      first_node_ep_ == SplineExtrapolation::linear ||
      last_node_ep_ == SplineExtrapolation::constant ||
      last_node_ep_ == SplineExtrapolation::linear;
    if (!needExtrapolationCoefficients)
        return;

    /* Beyond the spline nodes, we need to extrapolate using a * t + b.
     * Those coefficients are stored as
     * [
     *   D[b_first, p0], D[a_first, p0], D[b_last, p0], D[a_last, p0],
     *   D[b_first, p1], ... D[a_last, p{nplist-1}]
     * ]
     */
    coefficients_extrapolate_sensi.resize(4 * nplist, 0.0);

    realtype sm0;
    for (int ip = 0; ip < nplist; ip++) {
        realtype sp0 = dspline_valuesdp[spline_offset + ip];
        switch (first_node_ep_) {
            /* This whole switch-case-if-else-if-thing could be moved
             * outside the loop, I know. Yet, it's at most some thousand
             * if's done once in the program for saving many lines of code
             * and getting a much clearer code structure. */
            case SplineExtrapolation::constant:
                sm0 = 0;
                break;

            case SplineExtrapolation::linear:
                if (get_node_derivative_by_fd() &&
                    first_node_bc_ == SplineBoundaryCondition::given) {
                    sm0 =
                      (dspline_valuesdp[spline_offset + ip + nplist] - sp0) /
                      (nodes_[1] - nodes_[0]);

                } else if (get_node_derivative_by_fd() &&
                           first_node_bc_ == SplineBoundaryCondition::natural) {
                    throw AmiException(
                      "Natural boundary condition for "
                      "Hermite splines with linear extrapolation is "
                      "not yet implemented.");

                } else if (!get_node_derivative_by_fd() &&
                           first_node_bc_ == SplineBoundaryCondition::given) {
                    // sm0 = dspline_slopesdp[spline_offset + ip];
                    throw AmiException(
                      "Natural boundary condition for "
                      "Hermite splines with linear extrapolation is "
                      "not yet implemented.");

                } else if (!get_node_derivative_by_fd() &&
                           first_node_bc_ == SplineBoundaryCondition::natural) {
                    throw AmiException("Not implemented: sm0 is not set");

                } else {
                    throw AmiException(
                      "Some weird combination of spline boundary "
                      "condition, extrapolation and finite differences was "
                      "passed which should not be allowed.");
                }
                break;

            default:
                /* We don't need specific coefficients in the cases of:
                 * noExtrapolation, polynomial, periodic
                 *
                 * TODO: But we need to set sm0 even here,
                 *       otherwise its value is undefined!!!
                 *       Until corrected, raise an exception
                 */
                throw AmiException("Spline extrapolation sensitivity "
                                   "computation not supported yet "
                                   "for the following cases: noExtrapolation, "
                                   "polynomial, periodic");
        }
        /* Write them to the vector */
        coefficients_extrapolate_sensi[4 * ip] = sp0 - sm0 * nodes_[0];
        coefficients_extrapolate_sensi[4 * ip + 1] = sm0;
    }

    realtype sm_end;
    for (int ip = 0; ip < nplist; ip++) {
        realtype sp_end =
          dspline_valuesdp[spline_offset + ip + (n_nodes() - 1) * nplist];
        switch (last_node_ep_) {
            /* This whole switch-case-if-else-if-thing could be moved
             * outside the loop, I know. Yet, it's at most some thousand
             * if's done once in the program for saving many lines of code
             * and getting a much clearer code structure. */
            case SplineExtrapolation::constant:
                sm_end = 0;
                break;

            case SplineExtrapolation::linear:
                if (get_node_derivative_by_fd() &&
                    last_node_bc_ == SplineBoundaryCondition::given) {
                    sm_end =
                      (sp_end - dspline_valuesdp[spline_offset + ip +
                                                 (n_nodes() - 2) * nplist]) /
                      (nodes_[n_nodes() - 1] - nodes_[n_nodes() - 2]);

                } else if (get_node_derivative_by_fd() &&
                           last_node_bc_ == SplineBoundaryCondition::natural) {
                    throw AmiException(
                      "Natural boundary condition for "
                      "Hermite splines with linear extrapolation is "
                      "not yet implemented.");

                } else if (!get_node_derivative_by_fd() &&
                           last_node_bc_ == SplineBoundaryCondition::given) {
                    sm_end = dspline_slopesdp[spline_offset + ip +
                                              (n_nodes() - 1) * nplist];

                } else if (!get_node_derivative_by_fd() &&
                           last_node_bc_ == SplineBoundaryCondition::natural) {
                    throw AmiException(
                      "Natural boundary condition for "
                      "Hermite splines with linear extrapolation is "
                      "not yet implemented.");

                } else {
                    throw AmiException(
                      "Some weird combination of spline boundary "
                      "condition, extrapolation and finite differecnces was "
                      "passed which should not be allowed.");
                }
                break;

            default:
                /* We don't need specific coefficients in the cases of:
                 * noExtrapolation, polynomial, periodic
                 *
                 * TODO: But we need to set sm_end even here,
                 *       otherwise its value is undefined!!!
                 *       Until corrected, raise an exception
                 */
                throw AmiException("Spline extrapolation sensitivity "
                                   "computation not supported yet "
                                   "for the following cases: noExtrapolation, "
                                   "polynomial, periodic");
        }
        /* Write them to the vector */
        coefficients_extrapolate_sensi[4 * ip + 2] =
          sp_end - sm_end * nodes_[n_nodes() - 1];
        coefficients_extrapolate_sensi[4 * ip + 3] = sm_end;
    }
}

void
HermiteSpline::get_coeffs_sensi_lowlevel(int ip,
                                         int i_node,
                                         int nplist,
                                         int n_spline_coefficients,
                                         int spline_offset,
                                         realtype len,
                                         realtype len_m,
                                         realtype len_p,
                                         gsl::span<realtype> dnodesdp,
                                         gsl::span<realtype> dslopesdp,
                                         gsl::span<realtype> coeffs)
{
    /* We're using the short hand notation for node values and slopes from
     * compute_coefficients_sensi() here. See this function for documentation.
     */
    int node_offset = spline_offset + ip;
    int last = n_nodes() - 1;
    double spk = dnodesdp[node_offset + i_node * nplist];
    double spk1 = dnodesdp[node_offset + (i_node + 1) * nplist];
    double smk;
    double smk1;

    /* Get sensitivities of slopes. Depending on finite differences are used
     * or not, this may be a bit cumbersome now... */
    if (get_node_derivative_by_fd()) {
        if (i_node == 0) {
            /* Are we at the fist node? What's the boundary condition? */
            if (first_node_bc_ == SplineBoundaryCondition::zeroDerivative) {
                smk = 0;
            } else if (first_node_bc_ == SplineBoundaryCondition::given) {
                smk = (spk1 - spk) / len;
            } else if (first_node_bc_ == SplineBoundaryCondition::natural) {
                throw AmiException("Natural boundary condition for Hermite "
                                   "splines is not yet implemented.");
            } else if (first_node_bc_ == SplineBoundaryCondition::periodic) {
                smk = (spk1 - dnodesdp[node_offset + (last - 1) * nplist]) /
                      (len + nodes_[last] - nodes_[last - 1]);
            } else {
                /* must be SplineBoundaryCondition::naturalZeroDerivative*/
                throw AmiException(
                  "Natural boundary condition with zero "
                  "derivative is prohibited for Hermite splines.");
            }
            smk1 =
              (dnodesdp[node_offset + (i_node + 2) * nplist] - spk) / len_p;

        } else if (i_node == n_nodes() - 2) {
            /* Are we at the last node? What's the boundary condition? */
            smk =
              (spk1 - dnodesdp[node_offset + (i_node - 1) * nplist]) / len_m;
            if (last_node_bc_ == SplineBoundaryCondition::zeroDerivative) {
                smk1 = 0;
            } else if (last_node_bc_ == SplineBoundaryCondition::given) {
                smk1 = (spk1 - spk) / len;
            } else if (last_node_bc_ == SplineBoundaryCondition::natural) {
                throw AmiException("Natural boundary condition for Hermite "
                                   "splines is not yet implemented.");
            } else if (last_node_bc_ == SplineBoundaryCondition::periodic) {
                smk1 = (dnodesdp[node_offset + nplist] - spk) /
                       (len + nodes_[1] - nodes_[0]);
            } else {
                // must be SplineBoundaryCondition::naturalZeroDerivative
                throw AmiException(
                  "Natural boundary condition with zero "
                  "derivative is prohibited for Hermite splines.");
            }

        } else {
            /* We're somewhere in between. That's fine. */
            smk =
              (spk1 - dnodesdp[node_offset + (i_node - 1) * nplist]) / len_m;
            smk1 =
              (dnodesdp[node_offset + (i_node + 2) * nplist] - spk) / len_p;
        }
    } else {
        /* The slopes are explicitly given, easiest case... */
        smk = dslopesdp[node_offset + i_node * nplist];
        smk1 = dslopesdp[node_offset + (i_node + 1) * nplist];

        /* For the nodes at the boundary, we have to take care of the bc */
        if (i_node == 0 &&
            first_node_bc_ == SplineBoundaryCondition::zeroDerivative)
            smk = 0;
        if (i_node == n_nodes() - 2 &&
            last_node_bc_ == SplineBoundaryCondition::zeroDerivative)
            smk1 = 0;
    }

    /* Compute the actual coefficients */
    coeffs[ip * n_spline_coefficients + 4 * i_node] = spk;
    coeffs[ip * n_spline_coefficients + 4 * i_node + 1] = len * smk;
    coeffs[ip * n_spline_coefficients + 4 * i_node + 2] =
      3 * (spk1 - spk) - len * (2 * smk + smk1);
    coeffs[ip * n_spline_coefficients + 4 * i_node + 3] =
      2 * (spk - spk1) + len * (smk + smk1);
}

void
HermiteSpline::compute_final_value()
{
    /* We need to compute the final value of the spline, depending on its
     * boundary condition and the extrapolation option. */
    realtype finalValue;
    if ((last_node_ep_ == SplineExtrapolation::constant) ||
        (last_node_bc_ == SplineBoundaryCondition::zeroDerivative &&
         last_node_ep_ == SplineExtrapolation::linear)) {
        finalValue = coefficients_extrapolate[2];
    } else if (last_node_ep_ == SplineExtrapolation::linear) {
        if (coefficients_extrapolate[2] < 0) {
            finalValue = -INFINITY;
        } else if (coefficients_extrapolate[2] > 0) {
            finalValue = INFINITY;
        } else {
            finalValue = coefficients_extrapolate[2];
        }
    } else if (last_node_ep_ == SplineExtrapolation::polynomial) {
        int last = 4 * (n_nodes() - 1) - 1;
        if (coefficients[last] < 0) {
            finalValue = -INFINITY;
        } else if (coefficients[last] > 0) {
            finalValue = INFINITY;
        } else if (coefficients[last - 1] < 0) {
            finalValue = -INFINITY;
        } else if (coefficients[last - 1] > 0) {
            finalValue = INFINITY;
        } else if (coefficients[last - 2] < 0) {
            finalValue = -INFINITY;
        } else if (coefficients[last - 2] > 0) {
            finalValue = INFINITY;
        } else {
            finalValue = coefficients[last - 3];
        }
    } else {
        /* Periodic: will not yield a steady state */
        finalValue = NAN;
    }
    set_final_value(finalValue);
}

void
HermiteSpline::compute_final_sensitivity(
  int nplist,
  int /*spline_offset*/,
  gsl::span<realtype> /*dspline_valuesdp*/,
  gsl::span<realtype> /*dspline_slopesdp*/)
{
    /* We need to compute the final value of the spline, depending on its
     * boundary condition and the extrapolation option. */
    std::vector<realtype> finalSensitivity(nplist, 0);
    if ((last_node_ep_ == SplineExtrapolation::constant) ||
        (last_node_bc_ == SplineBoundaryCondition::zeroDerivative &&
         last_node_ep_ == SplineExtrapolation::linear)) {
        for (int ip = 0; ip < nplist; ip++)
            finalSensitivity[ip] = coefficients_extrapolate_sensi[4 * ip + 2];
    } else if (last_node_ep_ == SplineExtrapolation::linear) {
        /* If steady state is infinity, sensitivity must be 0 */
    } else if (last_node_ep_ == SplineExtrapolation::polynomial) {
        /* Yes, that's not correct. But I don't see any good reason for
         * implementing a case, which anybody with more than a dead fish
         * between the ears will never use. */
        std::fill(finalSensitivity.begin(), finalSensitivity.end(), NAN);
    } else {
        /* Periodic: will not yield a steady state */
        std::fill(finalSensitivity.begin(), finalSensitivity.end(), NAN);
    }
    set_final_sensitivity(finalSensitivity);
}

realtype
HermiteSpline::get_value(const double t) const
{
    /* Is this a steady state computation? */
    if (std::isinf(t))
        return get_final_value();

    /* Compute the spline value */
    int i_node;
    realtype len;

    /* Are we past the last node? Extrapolate! */
    if (t > nodes_[n_nodes() - 1]) {
        switch (last_node_ep_) {
            case SplineExtrapolation::noExtrapolation:
                throw AmiException(
                  "Trying to evaluate spline after last "
                  "spline node, but spline has been specified not to allow "
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
                return evaluate_polynomial(
                  (t - nodes_[i_node]) / len,
                  gsl::make_span(coefficients).subspan(i_node * 4));

            case SplineExtrapolation::periodic:
                len = nodes_[n_nodes() - 1] - nodes_[0];
                return get_value(nodes_[0] + std::fmod(t - nodes_[0], len));
            default:
                throw AmiException("Unsupported SplineExtrapolation type");
        }
    }

    /* Are we before the first node? Extrapolate! */
    if (t < nodes_[0]) {
        switch (first_node_ep_) {
            case SplineExtrapolation::noExtrapolation:
                throw AmiException(
                  "Trying to evaluate spline before first "
                  "spline node, but spline has been specified not to allow "
                  "extrapolation.");

            case SplineExtrapolation::constant:
                return coefficients_extrapolate[0];

            case SplineExtrapolation::linear:
                return coefficients_extrapolate[0] +
                       t * coefficients_extrapolate[1];

            case SplineExtrapolation::polynomial:
                /* Evaluate last interpolation polynomial */
                len = nodes_[1] - nodes_[0];
                return evaluate_polynomial((t - nodes_[0]) / len, coefficients);

            case SplineExtrapolation::periodic:
                len = nodes_[n_nodes() - 1] - nodes_[0];
                return get_value(nodes_[n_nodes() - 1] +
                                 std::fmod(t - nodes_[0], len));
            default:
                throw AmiException("Unsupported SplineExtrapolation type");
        }
    }

    /* Get the spline interval which we need */
    if (get_equidistant_spacing()) {
        /* equidistant spacing: just compute the interval */
        len = nodes_[1] - nodes_[0];
        i_node = static_cast<int>(std::trunc((t - nodes_[0]) / len));
        i_node = std::min(i_node, n_nodes() - 2);
    } else {
        /* no equidistant spacing: we need to iterate */
        i_node = 0;
        while (nodes_[i_node + 1] < t) {
            i_node++;
        }
        len = nodes_[i_node + 1] - nodes_[i_node];
    }

    /* Evaluate the interpolation polynomial */
    return evaluate_polynomial(
      (t - nodes_[i_node]) / len,
      gsl::make_span(coefficients).subspan(i_node * 4));
}

realtype
HermiteSpline::get_sensitivity(const double t, const int ip)
{
    /* Is this a steady state computation? */
    if (std::isinf(t))
        return get_final_sensitivity(ip);

    /* Compute the parametric derivative of the spline value */
    int i_node;
    realtype len;

    if (t > nodes_[n_nodes() - 1]) {
        /* Are we past the last node? Extrapolate! */
        switch (last_node_ep_) {
            case SplineExtrapolation::noExtrapolation:
                throw AmiException(
                  "Trying to evaluate spline sensitivity "
                  "after last spline node, but spline has been specified "
                  "to not allow extrapolation.");

            case SplineExtrapolation::constant:
                return coefficients_extrapolate_sensi[4 * ip + 2];

            case SplineExtrapolation::linear:
                return coefficients_extrapolate_sensi[4 * ip + 2] +
                       t * coefficients_extrapolate_sensi[4 * ip + 3];

            case SplineExtrapolation::polynomial:
                /* Evaluate last interpolation polynomial */
                i_node = n_nodes() - 2;
                len = nodes_[i_node + 1] - nodes_[i_node];
                return evaluate_polynomial(
                  (t - nodes_[i_node]) / len,
                  gsl::make_span(coefficients_sensi)
                    .subspan(ip * (n_nodes() - 1) * 4 + i_node * 4));

            case SplineExtrapolation::periodic:
                len = nodes_[n_nodes() - 1] - nodes_[0];
                return get_sensitivity(
                  nodes_[0] + std::fmod(t - nodes_[0], len), ip);
            default:
                throw AmiException("Unsupported SplineExtrapolation type");
        }
    }

    if (t < nodes_[0]) {
        /* Are we before the first node? Extrapolate! */
        switch (first_node_ep_) {
            case SplineExtrapolation::noExtrapolation:
                throw AmiException(
                  "Trying to evaluate spline before first "
                  "spline node, but spline has been specified to not allow "
                  "extrapolation.");

            case SplineExtrapolation::constant:
                return coefficients_extrapolate_sensi[4 * ip + 0];

            case SplineExtrapolation::linear:
                return coefficients_extrapolate_sensi[4 * ip + 0] +
                       t * coefficients_extrapolate_sensi[4 * ip + 1];

            case SplineExtrapolation::polynomial:
                /* Evaluate last interpolation polynomial */
                len = nodes_[1] - nodes_[0];
                return evaluate_polynomial(
                  (t - nodes_[0]) / len,
                  gsl::make_span(coefficients_sensi)
                    .subspan(ip * (n_nodes() - 1) * 4));

            case SplineExtrapolation::periodic:
                len = nodes_[n_nodes() - 1] - nodes_[0];
                return get_sensitivity(
                  nodes_[n_nodes() - 1] + std::fmod(t - nodes_[0], len), ip);
            default:
                throw AmiException("Unsupported SplineExtrapolation type");
        }
    }

    /* Get the spline interval which we need */
    if (get_equidistant_spacing()) {
        /* equidistant spacing: just compute the interval */
        len = nodes_[1] - nodes_[0];
        i_node = static_cast<int>(std::trunc((t - nodes_[0]) / len));
        i_node = std::min(i_node, n_nodes() - 2);
    } else {
        /* no equidistant spacing: we need to iterate */
        i_node = 0;
        while (nodes_[i_node + 1] < t) {
            i_node++;
        }
        len = nodes_[i_node + 1] - nodes_[i_node];
    }

    /* Evaluate the interpolation polynomial */
    return evaluate_polynomial(
      (t - nodes_[i_node]) / len,
      gsl::make_span(coefficients_sensi)
        .subspan(ip * (n_nodes() - 1) * 4 + i_node * 4));
}

} // namespace amici
