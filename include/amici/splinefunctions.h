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
    explicit SplineFunction(const Model &model);
  
    ~SplineFunction(){};
    
    virtual void computeSpline() = 0;
    
    virtual double getValue(const double t) = 0;
    
    virtual double getParametricDerivative(const double t, const int ip) = 0;

    /**
     * @brief Accessor to equidistant_spacing_ member
     * @return equidistant_spacing flag
     */
    bool get_equidistant_spacing() {
        return equidistant_spacing_;
    }
    
    /**
     * @brief Accessor to node_derivative_by_FD_ member
     * @return node_derivative_by_FD flag
     */
    bool get_node_derivative_by_FD() {
        return node_derivative_by_FD_;
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
     * @brief Switch computing of derivatives at spline nodes by 
     * finite differences on or off 
     * @param node_derivative_by_FD flag for usage of finite differences
     */
    void set_node_derivative_by_FD(bool node_derivative_by_FD) {
        node_derivative_by_FD_ = node_derivative_by_FD;
    }

    /**
     * @brief Switch enforced positivity by logarithmic parametrization 
     * on or off 
     * @param logarithmic_paraterization flag for logarithmic parametrization
     */
    void set_logarithmic_paraterization(bool logarithmic_paraterization) {
        logarithmic_paraterization_ = logarithmic_paraterization;
    }    
    
  private:
    std::vector<realtype> nodes_;
    
    std::vector<realtype> node_values_;
    
    std::vector<realtype> node_values_derivative_;
    
    AmiVector coeffictions_;
    
    AmiVectorArray coeffictions_sensi_;
    
    bool node_derivative_by_FD_ = false;
    
    bool equidistant_spacing_ = false;
    
    bool logarithmic_paraterization_ = false;
    

}; // class SplineFunction



class HermiteSpline : splineFunction()

} // namespace amici

#endif /* AMICI_SPLINEFUNCTION_H */