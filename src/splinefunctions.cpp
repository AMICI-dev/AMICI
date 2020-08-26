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
    : nodes(nodes), node_values(node_values), 
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
        nodes.resize(n_nodes_);
        for (int i_node = 0; i_node < n_nodes(); i_node++)
            nodes[i_node] = node_start + i_node * node_step;
    }
}

HermiteSpline::HermiteSpline(std::vector<realtype> nodes,
                             std::vector<realtype> node_values,
                             std::vector<realtype> node_values_derivative,
                             SplineBoundaryCondition firstNodeDerivative,
                             SplineBoundaryCondition lastNodeDerivative,
                             bool node_derivative_by_FD,
                             bool equidistant_spacing,
                             bool logarithmic_paraterization) 
    : AbstractSpline(nodes, node_values, equidistant_spacing,
                     logarithmic_paraterization),
    node_values_derivative(std::move(node_values_derivative)),
    firstNodeDerivative(firstNodeDerivative), 
    lastNodeDerivative(lastNodeDerivative),
    node_derivative_by_FD_(node_derivative_by_FD) {
    
    /* If values of the derivative at the nodes are to be computed by finite 
     * differences, we have to fill up node_values_derivative_ */
    if (node_derivative_by_FD_) {
        node_values_derivative.resize(n_nodes(), 0.0);
        for (int i_node = 1; i_node < n_nodes() - 2; i_node++)
            node_values_derivative[i_node] = 
                (node_values[i_node + 1] - node_values[i_node - 1]) /
                (nodes[i_node + 1] - nodes[i_node - 1]);

        /* We have to take care of the first node (1-sided FD)... */
        switch (firstNodeDerivative) {
            case SplineBoundaryCondition::constant:
                node_values_derivative[0] = 0;
                break;

            case SplineBoundaryCondition::linearFinDiff:
                node_values_derivative[0] =
                    (node_values[1] - node_values[0]) / (nodes[1] - nodes[0]);
                break;
            
            case SplineBoundaryCondition::linearNatural:
                throw AmiException("Natural boundary condition for Hermite splines "
                                   "is not yet implemented.");
        }
        /* ...and the last node (1-sided FD). */
        switch (lastNodeDerivative) {
            case SplineBoundaryCondition::constant:
                node_values_derivative[n_nodes() - 1] = 0;
                break;

            case SplineBoundaryCondition::linearFinDiff:
                node_values_derivative[n_nodes() - 1] =
                    (node_values[n_nodes() - 1] - node_values[n_nodes() - 2]) / 
                    (nodes[n_nodes() - 1] - nodes[n_nodes() - 2]);
                break;
            
            case SplineBoundaryCondition::linearNatural:
                throw AmiException("Natural boundary condition for Hermite splines "
                                   "is not yet implemented.");
        }
    }
    
}

void HermiteSpline::computeCoefficients() {
    /* Allocate space for the coefficients for Horner's method. 
     * They are stored in the vector as 
     * [d_0, c_0, b_0, a_0, d_1, c_1, ... , b_{n_nodes-1}, a_{n_nodes-1}] */
    coefficients.resize(4 * (n_nodes() - 1), 0.0);
    /* Beyond the spline nodes, we need to extrapolate using a * t + b. 
     * Those coefficients are stored as [b_first, a_first, b_last, a_last] */
    coefficients_extrapolate.resize(4, 0.0);
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
        len = nodes[i_node + 1] - nodes[i_node];
        
        /* Coefficients for cubic Hermite polynomials */
        coefficients[4 * i_node] = node_values[i_node];
        coefficients[4 * i_node + 1] = len*node_values_derivative[i_node];
        coefficients[4 * i_node + 2] = - 3 * node_values[i_node]
            - 2 * len * node_values_derivative[i_node]
            + 3 * node_values[i_node + 1] 
            - len * node_values_derivative[i_node + 1];
        coefficients[4 * i_node + 3] = 2 * node_values[i_node]
            + len * node_values_derivative[i_node]
            - 2 * node_values[i_node + 1]
            + len * node_values_derivative[i_node + 1];
    }
    
    /* Coefficients for affine functions for extrapolation */
    coefficients_extrapolate[0] = node_values[0]
        - node_values[0] * node_values_derivative[0];
    coefficients_extrapolate[1] = node_values_derivative[0];
    coefficients_extrapolate[2] = node_values[n_nodes() - 1]
        - node_values[n_nodes() - 1] * node_values_derivative[n_nodes() - 1];
    coefficients_extrapolate[3] = node_values_derivative[n_nodes() - 1];
}

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
    /* Beyond the spline nodes, we need to extrapolate using a * t + b. 
     * Those coefficients are stored as 
     * [ D[b_first, p0], D[a_first, p0], D[b_last, p0], D[a_last, p0],
     *   D[b_first, p1], ... D[a_last, p{nplist-1}]
     * ] */
    coefficients_extrapolate_sensi.resize(4 * nplist, 0.0);
    realtype len, len_m, len_p;
 
    /* Parametric derivatives of splines are splines again.
     * We compute the coefficients for those polynomials now. */    
    for (int i_node = 0; i_node < n_nodes() - 2; i_node++) {
        /* Get the length of the interval. */
        len = nodes[i_node + 1] - nodes[i_node];
        len_m = (i_node > 0) ? nodes[i_node + 1] - nodes[i_node - 1] : len;
        len_p = (i_node < n_nodes() - 1) ? nodes[i_node + 2] - nodes[i_node] : len;

        /* As computing the coefficient is a mess, it's in another function */
        for (int ip = 0; ip < nplist; ip++)
            getCoeffsSensiLowlevel(ip, i_node, n_spline_coefficients, 
                                   spline_offset, len, len_m, len_p, 
                                   dspline_valuesdp, dspline_slopesdp, 
                                   coefficients_sensi.data(),
                                   coefficients_extrapolate_sensi.data());
    }

    /* We need the coefficients for extrapolating beyond the spline domain */
    for (int ip = 0; ip < nplist; ip++) {
        /* Before node[0] */
        realtype sp0 = dspline_valuesdp[n_nodes() * ip];
        realtype sm0;
        if (get_node_derivative_by_FD()) {
            if (firstNodeDerivative == SplineBoundaryCondition::constant) {
                sm0 = 0;
            } else if (firstNodeDerivative == SplineBoundaryCondition::linearFinDiff) {
                sm0 = (dspline_valuesdp[ip] - sp0) / (nodes[1] - nodes[0]);
            } else {
                throw AmiException("Natural boundary condition for Hermite splines "
                                   "is not yet implemented.");
            }
        } else {
            sm0 = dspline_slopesdp[n_nodes() * ip];
        }

        /* After node[n_nodes() - 1] */
        realtype sp_end = dspline_valuesdp[n_nodes() * (ip + 1) - 1];
        realtype sm_end;
        if (get_node_derivative_by_FD()) {
            if (firstNodeDerivative == SplineBoundaryCondition::constant) {
                sm_end = 0;
            } else if (firstNodeDerivative == SplineBoundaryCondition::linearFinDiff) {
                sm_end = (sp_end - dspline_valuesdp[n_nodes() * (ip + 1) - 1]) 
                    / (nodes[n_nodes() - 1] - nodes[n_nodes() - 2]);
            } else {
                throw AmiException("Natural boundary condition for Hermite splines "
                                   "is not yet implemented.");
            }
        } else {
            sm_end = dspline_slopesdp[n_nodes() * (ip + 1) - 1];
        }

        /* Write them to the vector */
        coefficients_extrapolate_sensi[4 * ip] = sp0 - sm0 * nodes[0];
        coefficients_extrapolate_sensi[4 * ip + 1] = sm0;
        coefficients_extrapolate_sensi[4 * ip + 2] = sp_end - sm_end * nodes[n_nodes() - 1];
        coefficients_extrapolate_sensi[4 * ip + 3] = sm_end;
    }
}

void HermiteSpline::getCoeffsSensiLowlevel(int ip, int i_node, int n_spline_coefficients, 
                                           int spline_offset, realtype len, realtype len_m, 
                                           realtype len_p, realtype *dnodesdp, 
                                           realtype *dslopesdp, realtype *coeffs, 
                                           realtype *coeffs_extrapol) {
    int node_offset = spline_offset + n_nodes() * ip;
    double spk = dnodesdp[node_offset + i_node];
    double spk1 = dnodesdp[node_offset + i_node + 1];
    double smk;
    double smk1; 

    /* Get sensitivities of slopes. Depending on finite differences are used 
     * or not, this may be a bit cumbersome now... */
    if (get_node_derivative_by_FD()) {
        if (i_node == 0) {
            /* Are we at the fist node? What's the boundary condition? */
            if (firstNodeDerivative == SplineBoundaryCondition::constant) {
                smk = 0;
            } else if (firstNodeDerivative == SplineBoundaryCondition::linearFinDiff) {
                smk = (spk1 - spk) / len_m;
            } else {
                throw AmiException("Natural boundary condition for Hermite splines "
                                   "is not yet implemented.");
            }
            smk1 = (dnodesdp[node_offset + i_node + 2] - spk) / len_p;
            
        } else if (i_node == n_nodes() - 1) {
            /* Are we at the last node? What's the boundary condition? */
            smk = (spk1 - dnodesdp[node_offset + i_node - 1]) / len_m;
            if (lastNodeDerivative == SplineBoundaryCondition::constant) {
                smk1 = 0;
            } else if (lastNodeDerivative == SplineBoundaryCondition::linearFinDiff) {
                smk1 = (spk1 - spk) / len_p;
            } else {
                throw AmiException("Natural boundary condition for Hermite splines "
                                   "is not yet implemented.");
            }

        } else {
            /* We're somewhere in between. That's fine. */
            smk = (spk1 - dnodesdp[node_offset + i_node - 1]) / len_m;
            smk1 = (dnodesdp[node_offset + i_node + 2] - spk) / len_p;
        }
    } else {
        /* The slopes are explicitly given, easiest case... */
        smk = dslopesdp[node_offset + i_node];
        smk1 = dslopesdp[node_offset + i_node + 1];
    }

    /* Compute the actual coefficients */
    coeffs[ip * n_spline_coefficients + 4 * i_node] = spk;
    coeffs[ip * n_spline_coefficients + 4 * i_node + 1] = len * smk;;
    coeffs[ip * n_spline_coefficients + 4 * i_node + 2] = 3 * (spk1  - spk) - len * (2 * smk + smk1);
    coeffs[ip * n_spline_coefficients + 4 * i_node + 3] = 2 * (spk - spk1) + len * (smk + smk1);
}

realtype HermiteSpline::getValue(const double t) {
    /* Compute the spline value */
    int i_node = 0;
    realtype len = nodes[1] - nodes[0];
    
    /* Are we past the last node? Extrapolate! */
    if (t > nodes[n_nodes() - 1])
        return coefficients_extrapolate[2] + t * coefficients_extrapolate[3];

    /* Are we before the first node? Extrapolate! */
    if (t < nodes[0])
        return coefficients_extrapolate[0] + t * coefficients_extrapolate[1];

    /* Get the spline interval which we need */
    if (get_equidistant_spacing()) {
        /* equidistant spacing: just compute the interval */
        i_node = std::floor((t - nodes[0]) / len);
    } else {
        /* no equidistant spacing: we need to iterate */
        while (nodes[i_node + 1] < t) {
            i_node++;
        }
        len = nodes[i_node + 1] - nodes[i_node];
    }

    /* Evaluate the interpolation polynomial */
    return evaluatePolynomial((t - nodes[i_node]) / len, 
                              &(coefficients[i_node * 4]));
}

realtype HermiteSpline::getSensitivity(const double t, const int ip) {
    /* Compute the parametric derivative of the spline value */
    if (t > nodes[n_nodes() - 1]) {
        /* Are we past the last node? Extrapolate! */
        return coefficients_extrapolate_sensi[4 * ip + 2]
            + t * coefficients_extrapolate_sensi[4 * ip + 3];

    } else if (t < nodes[0]) {
        /* Are we before the first node? Extrapolate! */
        return coefficients_extrapolate_sensi[4 * ip + 0]
            + t * coefficients_extrapolate_sensi[4 * ip + 1];

    } else {
        /* Get the spline interval which we need */
        realtype len;
        int i_node = 0;
        if (get_equidistant_spacing()) {
            /* equidistant spacing: just compute the interval */
            len = nodes[1] - nodes[0];
            i_node = std::floor((t - nodes[0]) / len);
        } else {
            /* no equidistant spacing: we need to iterate */
            while (nodes[i_node + 1] < t)
                i_node++;

            len = nodes[i_node + 1] - nodes[i_node];
        }

        return evaluatePolynomial((t - nodes[i_node]) / len, 
                                  &coefficients_sensi[ip * n_nodes() * 4 + i_node * 4]);
    }
    
}

} // namespace amici
