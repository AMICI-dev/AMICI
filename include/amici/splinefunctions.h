#ifndef AMICI_SPLINEFUNCTION_H
#define AMICI_SPLINEFUNCTION_H

#include "amici/defines.h"
#include "amici/vector.h"

namespace amici {
    
class Model;
/**
 * @brief The spline class is an AMICI-own implementation. Instances of this 
 * class are created upon solver setup and the needed splines are set up 
 * (e.g., interpolation of the nodes is performed). 
 * Upon call to a spline fuction, only the evaluation of the spline polynomial
 * is carried out.
 */
class SplineFunction {
  public:
    /**
     * @brief constructor
     * @param model Model instance
     */
    explicit SplineFunction(std::vector<realtype> nodes, 
                            std::vector<realtype> node_values,
                            bool equidistant_spacing,
                            bool logarithmic_paraterization);
  
    ~SplineFunction(){};
    
    virtual void computeCoefficients() = 0;
    
    virtual void computeCoefficientsSensi() = 0;
    
    virtual double getValue(const double t) = 0;
    
    virtual double getSensitivity(const double t, const int ip) = 0;

    /**
     * @brief Accessor to equidistant_spacing_ member
     * @return equidistant_spacing flag
     */
    bool get_equidistant_spacing() {
        return equidistant_spacing_;
    }

    /**
     * @brief Accessor to logarithmic_paraterization_ member
     * @return logarithmic_paraterization flag
     */
    bool get_logarithmic_paraterization() {
        return logarithmic_paraterization_;
    }

  protected: 
  /*
   * In order to have the main data members private, we need protected 
   * accessor macros.
   * */
  
    /**
     * @brief Switch equisitant spacing of spline nodes on or off 
     * @param equidistant_spacing flag for equidistancy of spline nodes
     */
    void set_equidistant_spacing(bool equidistant_spacing) {
        equidistant_spacing_ = equidistant_spacing;
    }

    /**
     * @brief Switch enforced positivity by logarithmic parametrization 
     * on or off 
     * @param logarithmic_paraterization flag for logarithmic parametrization
     */
    void set_logarithmic_paraterization(bool logarithmic_paraterization) {
        logarithmic_paraterization_ = logarithmic_paraterization;
    }

    std::vector<realtype> nodes;

    std::vector<realtype> node_values;

    std::vector<realtype> coefficients;
    
    std::vector<realtype> coefficients_extrapolate;

    std::vector<realtype> coefficients_sensi;
    
    std::vector<realtype> coefficients_extrapolate_sensi;

    const int n_nodes() { return n_nodes; }
    
  private:
    bool equidistant_spacing_ = false;
    
    bool logarithmic_paraterization_ = false;
    
    int n_nodes;
    
}; // class SplineFunction



class HermiteSpline : splineFunction() {
  public:
    HermiteSpline() = default;
      
    HermiteSpline(std::vector<realtype> nodes,
                  std::vector<realtype> node_values,
                  std::vector<realtype> node_values_derivative,
                  SplineBoundaryCondition firstNodeDerivative,
                  SplineBoundaryCondition lastNodeDerivative,
                  bool node_derivative_by_FD,
                  bool equidistant_spacing,
                  bool logarithmic_paraterization);

    ~HermiteSpline(){};
    
    void computeCoefficients() override;
    
    void computeCoefficientsSensi(Model *model) override;
    
    double getValue(const double t) override;
    
    double getSensitivity(const double t, const int ip) override;
    
    /**
     * @brief Accessor to node_derivative_by_FD_ member
     * @return node_derivative_by_FD flag
     */
    bool get_node_derivative_by_FD() {
        return node_derivative_by_FD_;
    }

  private:
    void getCoeffsSensiLowlevel(int ip, int i_node, int offset, 
                                realtype len, realtype, len_m, realtype len_p,
                                realtype *dnodesdp, realtype *dslopesdp,
                                realtype *coeffs, realtype *coeffs_extrapol);

    std::vector<realtype> node_values_derivative;
      
    bool node_derivative_by_FD_ = false;
    
    SplineBoundaryCondition firstNodeDerivative = SplineBoundaryCondition::linearFinDiff;
    
    SplineBoundaryCondition lastNodeDerivative = SplineBoundaryCondition::constant;
    
}; // class HermiteSpline

} // namespace amici

#endif /* AMICI_SPLINEFUNCTION_H */