#include "amici/splinefunctions.h"
#include "amici/defines.h"
#include "amici/vector.h"
#include "amici/amici.h"
#include "amici/exception.h"

namespace amici {
    
SplineFunction::SplineFunction(std::vector<realtype> nodes, 
                               std::vector<realtype> node_values,
                               bool equidistant_spacing,
                               bool logarithmic_paraterization) 
    : nodes(nodes), node_values(node_values), 
    equidistant_spacing_(equidistant_spacing),
    logarithmic_paraterization_(logarithmic_paraterization) {
    
    /* we want to set the number of nodes and compute node locations, 
     * in case we have equidistant spacing */
    n_nodes = node_values.size();
    if (equidistant_spacing_) {
        if (nodes.size() != 2)
            throw AmiException("Splines with equidistant spacing need a nodes "
                               "vector with two elements (first/last node).");
        realtype node_start =  nodes[0];
        realtype node_step = (nodes[1] - nodes[0]) / (n_nodes - 1);
        nodes.resize(n_nodes);
        for (i_node = 0; i_node < n_nodes; i_node++)
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
        for (i_node = 1; i_node < n_nodes() - 2; i_node++)
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

} // namespace amici
