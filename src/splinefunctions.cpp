#include "amici/splinefunctions.h"
#include "amici/defines.h"
#include "amici/vector.h"
#include "amici/amici.h"
#include "amici/exception.h"

namespace amici {

static realtype evaluatePolynomial(realtype x, realtype coeff[4]) {
    return coeff[0] + x * (coeff[1] + x * (coeff[2] + x * coeff[3]));
}
    
SplineFunction::SplineFunction(std::vector<realtype> nodes, 
                               std::vector<realtype> node_values,
                               bool equidistant_spacing,
                               bool logarithmic_paraterization,
                               SensitivityMethod sensi) 
    : nodes(nodes), node_values(node_values), 
    equidistant_spacing_(equidistant_spacing),
    logarithmic_paraterization_(logarithmic_paraterization) {
    
    /* we want to set the number of nodes */
    n_nodes = node_values.size();
    
    /* In case we have equidistant spacing, compute node locations */
    if (equidistant_spacing_) {
        if (nodes.size() != 2)
            throw AmiException("Splines with equidistant spacing need a nodes "
                               "vector with two elements (first/last node).");
        realtype node_start =  nodes[0];
        realtype node_step = (nodes[1] - nodes[0]) / (n_nodes - 1);
        nodes.resize(n_nodes);
        for (int i_node = 0; i_node < n_nodes; i_node++)
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
    : splineFunction(nodes, node_values, equidistant_spacing,
                     logarithmic_paraterization),
    node_values_derivative(node_values_derivative),
    firstNodeDerivative(firstNodeDerivative), 
    lastNodeDerivative(lastNodeDerivative),
    node_derivative_by_FD_(node_derivative_by_FD) {
    
    /* If values of the derivative at the nodes are to be computed by finite 
     * differences, we have to fill up node_values_derivative_ */
    if (node_derivative_by_FD_) {
        node_values_derivative.resize(n_nodes, 0.0);
        for (int i_node = 1; i_node < n_nodes() - 2; i_node++)
            node_values_derivative[i_node] = 
                (node_values[i_node + 1] - node_values[i_node - 1]) /
                (nodes[i_node + 1] - nodes[i_node - 1]);
    }
    
    /* We have to take care of the first node... */
    switch (firstNodeDerivative) {
        case SplineBoundaryCondition.constant:
            node_values_derivative[0] = 0;
            break;

        case SplineBoundaryCondition.linearFinDiff:
            node_values_derivative[0] =
                (node_values[1] - node_values[0]) / (nodes[1] - nodes[0]);
            break;
        
        case SplineBoundaryCondition.linearNatural:
            throw AmiException("Natural boundary condition for Hermite splines "
                               "is not yet implemented.");
    }
    /* ...and the last node. */
    switch (lastNodeDerivative) {
        case SplineBoundaryCondition.constant:
            node_values_derivative[n_nodes() - 1] = 0;
            break;

        case SplineBoundaryCondition.linearFinDiff:
            node_values_derivative[n_nodes() - 1] =
                (node_values[n_nodes() - 1] - node_values[n_nodes() - 2]) / 
                (nodes[n_nodes() - 1] - nodes[n_nodes() - 2]);
            break;
        
        case SplineBoundaryCondition.linearNatural:
            throw AmiException("Natural boundary condition for Hermite splines "
                               "is not yet implemented.");
    }
    
}

HermiteSpline::computeCoefficients() {
    /* allocate space for the coefficients */
    coefficients.resize(4 * (n_nodes() - 1), 0.0);
    
    /* Compute the coefficients of the spline polynomials when 
     * using Horner's method for each interval:
     * Horner's method: 
     * spline(t) = a * t**3 + b * t**2 + c * t + d
     *           = d + t * (c + t * (b + t * a))
     * with coefficients[4 * i_node + (0, 1, 2, 3)] = (d, c, b, a)
     * */
    realtype len;
    realtype tmp_coeff[4];
 
    for (int i_node = 0; i_node < n_nodes() - 1; i_node++) {
        len = nodes[i_node + 1] - nodes[i_node];
        coefficients[4 * i_node] = node_values[i_node];
        coefficients[4 * i_node + 1] = len*node_values_derivative[i_node];
        coefficients[4 * i_node + 2] = - 3 * node_values[i_node]
            - 2 * len * node_values_derivative[i_node]
            + 3 * node_values[i_node + 1] 
            - len * node_values_derivative[i_node + 1];
        coefficients[4 * i_node + 3] = 2 * node_values[i_node]
            + len * node_values_derivative[i_node]
            - 2 * node_values[i_node + 1
            + len * node_values_derivative[i_node + 1];
    }
}

HermiteSpline::getValue(const double t) {
    /* Compute the spline value */
    if (t > nodes[n_nodes() - 1]) {
        /* Are we past the last node? Extrapolate! */
        if (firstNodeDerivative == SplineBoundaryCondition::constant) {
            return node_values[n_nodes() - 1];
        } else {
            return node_values[n_nodes() - 1] + 
                (t - node_values[n_nodes() - 1]) * 
                node_values_derivative[n_nodes() - 1];
        }

    } else if (t < nodes[0]) {
        /* Are we before the first node? Extrapolate! */
        if (firstNodeDerivative == SplineBoundaryCondition::constant) {
            return node_values[0];
        } else {
            return node_values[0] + (t - node_values[0]) * 
                node_values_derivative[0];
        }

    } else {
        /* Get the spline interval which we need */
        if (get_equidistant_spacing) {
            /* equidistant spacing: just compute the interval */
            realtype len = nodes[1] - nodes[0];
            int i_node = std::floor((t - nodes[0]) / len);
        } else {
            /* no equidistant spacing: we need to iterate */
            int i_node = 0;
            while (nodes[i_node + 1] < t)
                i_node++;

            realtype len = nodes[i_node + 1] - nodes[i_node];
        }

        return evaluatePolynomial((t - nodes[i_node]) / len, 
                                  &coefficients[i_node * 4]);
    }
    
}

} // namespace amici
