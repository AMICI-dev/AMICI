#ifndef AMICI_SPLINEFUNCTIONS_H
#define AMICI_SPLINEFUNCTIONS_H

#include "amici/defines.h"
#include <vector>

namespace amici {
    
class Model;
/**
 * @brief The spline class is an AMICI-own implementation. Instances of this 
 * class are created upon solver setup and the needed splines are set up 
 * (e.g., interpolation of the nodes is performed). 
 * Upon call to a spline fuction, only the evaluation of the spline polynomial
 * is carried out.
 */
class AbstractSpline {
  public:
    /** default constructor */
    AbstractSpline() = default;

    /**
     * @brief constructor
     * @param model Model instance
     */
    AbstractSpline(std::vector<realtype> nodes, 
                   std::vector<realtype> node_values,
                   bool equidistant_spacing,
                   bool logarithmic_paraterization);
  
    ~AbstractSpline(){};
    
    virtual void computeCoefficients() = 0;
    
    virtual void computeCoefficientsSensi(int nplist, int spline_offset, 
                                          realtype *dnodesdp, 
                                          realtype *dslopesdp) = 0;
    
    virtual double getValue(const double t) = 0;
    
    virtual double getSensitivity(const double t, const int ip) = 0;

    virtual bool get_node_derivative_by_FD() = 0;

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

    const int n_nodes() { return n_nodes_; }

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
    
  private:
    bool equidistant_spacing_ = false;
    
    bool logarithmic_paraterization_ = false;
    
    int n_nodes_;
    
}; // class SplineFunction



class HermiteSpline : public AbstractSpline {
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
    
    void computeCoefficientsSensi(int nplist, int spline_offset, 
                                  realtype *dnodesdp,
                                  realtype *dslopesdp) override;
    
    double getValue(const double t) override;
    
    double getSensitivity(const double t, const int ip) override;
    
    /**
     * @brief Accessor to node_derivative_by_FD_ member
     * @return node_derivative_by_FD flag
     */
    bool get_node_derivative_by_FD() override {
        return node_derivative_by_FD_;
    }

  private:
    void getCoeffsSensiLowlevel(int ip, int i_node, int n_spline_coefficients, 
                                int spline_offset, realtype len, realtype len_m, 
                                realtype len_p, realtype *dnodesdp, realtype *dslopesdp,
                                realtype *coeffs, realtype *coeffs_extrapol);

    std::vector<realtype> node_values_derivative;

    SplineBoundaryCondition firstNodeDerivative = SplineBoundaryCondition::linearFinDiff;
    
    SplineBoundaryCondition lastNodeDerivative = SplineBoundaryCondition::constant;

    bool node_derivative_by_FD_ = false;
    
}; // class HermiteSpline

} // namespace amici

#endif /* AMICI_SPLINEFUNCTION_H */