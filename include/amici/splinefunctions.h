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
                   bool logarithmic_parametrization);

    virtual ~AbstractSpline(){};

    virtual void computeCoefficients() = 0;

    virtual void computeCoefficientsSensi(int nplist, int spline_offset,
                                          realtype *dnodesdp,
                                          realtype *dslopesdp) = 0;

    virtual realtype getValue(const realtype t) = 0;

    virtual realtype getSensitivity(const realtype t, const int ip) = 0;

    virtual bool get_node_derivative_by_FD() = 0;

    /**
     * @brief Accessor to equidistant_spacing_ member
     * @return equidistant_spacing flag
     */
    bool get_equidistant_spacing();

    /**
     * @brief Accessor to logarithmic_parametrization_ member
     * @return logarithmic_parametrization flag
     */
    bool get_logarithmic_parametrization();

    const int n_nodes() { return n_nodes_; }

  protected:
    std::vector<realtype> nodes_;

    std::vector<realtype> node_values_;

    std::vector<realtype> coefficients;

    std::vector<realtype> coefficients_extrapolate;

    std::vector<realtype> coefficients_sensi;

    std::vector<realtype> coefficients_extrapolate_sensi;

    virtual void computeFinalValue() = 0;

    virtual void computeFinalSensitivity(int nplist, int spline_offset,
                                         realtype *dspline_valuesdp,
                                         realtype *dspline_slopesdp) = 0;

    realtype getFinalValue();

    void setFinalValue(realtype finalValue);

    realtype getFinalSensitivity(const int ip);

    void setFinalSensitivity(std::vector<realtype> const finalSensitivity);

  /*
   * In order to have the main data members private, we need protected
   * accessor macros.
   * */

    /**
     * @brief Switch equisitant spacing of spline nodes on or off
     * @param equidistant_spacing flag for equidistancy of spline nodes
     */
    void set_equidistant_spacing(bool equidistant_spacing);

    /**
     * @brief Switch enforced positivity by logarithmic parametrization
     * on or off
     * @param logarithmic_parametrization flag for logarithmic parametrization
     */
    void set_logarithmic_parametrization(bool logarithmic_parametrization);

  private:
    realtype finalValue_;

    std::vector<realtype> finalSensitivity_;

    bool equidistant_spacing_ = false;

    bool logarithmic_parametrization_ = false;

    int n_nodes_;

}; // class SplineFunction



class HermiteSpline : public AbstractSpline {
  public:
    HermiteSpline() = default;

    HermiteSpline(std::vector<realtype> nodes,
                  std::vector<realtype> node_values,
                  std::vector<realtype> node_values_derivative,
                  SplineBoundaryCondition firstNodeBC,
                  SplineBoundaryCondition lastNodeBC,
                  SplineExtrapolation firstNodeExtrapol,
                  SplineExtrapolation lastNodeExtrapol,
                  bool node_derivative_by_FD,
                  bool equidistant_spacing,
                  bool logarithmic_parametrization);

    ~HermiteSpline(){};

    void computeCoefficients() override;

    void computeCoefficientsSensi(int nplist, int spline_offset,
                                  realtype *dnodesdp,
                                  realtype *dslopesdp) override;

    realtype getValue(const double t) override;

    realtype getSensitivity(const double t, const int ip) override;

    /**
     * @brief Accessor to node_derivative_by_FD_ member
     * @return node_derivative_by_FD flag
     */
    bool get_node_derivative_by_FD() override {
        return node_derivative_by_FD_;
    }

  private:
    void getCoeffsSensiLowlevel(int ip, int i_node, int nplist, int n_spline_coefficients,
                                int spline_offset, realtype len, realtype len_m,
                                realtype len_p, realtype *dnodesdp, realtype *dslopesdp,
                                realtype *coeffs, realtype *coeffs_extrapol);

    void handleInnerDerviatives();

    void handleBoundaryConditions();

    void computeCoefficientsExtrapolation();

    void computeCoefficientsExtrapolationSensi(int nplist, int spline_offset,
                                               realtype *dspline_valuesdp,
                                               realtype *dspline_slopesdp);

    void computeFinalValue() override;

    void computeFinalSensitivity(int nplist, int spline_offset,
                                 realtype *dspline_valuesdp,
                                 realtype *dspline_slopesdp) override;

    std::vector<realtype> node_values_derivative_;

    SplineBoundaryCondition firstNodeBC_ = SplineBoundaryCondition::given;

    SplineBoundaryCondition lastNodeBC_ = SplineBoundaryCondition::given;

    SplineExtrapolation firstNodeEP_ = SplineExtrapolation::linear;

    SplineExtrapolation lastNodeEP_ = SplineExtrapolation::linear;

    bool node_derivative_by_FD_ = false;

}; // class HermiteSpline

} // namespace amici

#endif /* AMICI_SPLINEFUNCTION_H */
