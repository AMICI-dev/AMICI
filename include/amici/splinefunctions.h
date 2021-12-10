#ifndef AMICI_SPLINEFUNCTIONS_H
#define AMICI_SPLINEFUNCTIONS_H

#include "amici/defines.h"

#include <vector>

#include <gsl/gsl-lite.hpp>

namespace amici {

class Model;
/**
 * @brief AMICI spline base class.
 *
 * Instances of this class are created upon solver setup and the needed splines
 * are set up (e.g., interpolation of the nodes is performed).
 * Upon call to a spline fuction, only the evaluation of the spline polynomial
 * is carried out.
 */
class AbstractSpline
{
  public:
    /** default constructor */
    AbstractSpline() = default;

    /**
     * @brief AbstractSpline TODO
     * @param nodes
     * @param node_values
     * @param equidistant_spacing
     * @param logarithmic_parametrization
     */
    AbstractSpline(std::vector<realtype> nodes,
                   std::vector<realtype> node_values,
                   bool equidistant_spacing,
                   bool logarithmic_parametrization);

    virtual ~AbstractSpline() = default;

    virtual void compute_coefficients() = 0;

    virtual void compute_coefficients_sensi(int nplist,
                                            int spline_offset,
                                            gsl::span<realtype> dnodesdp,
                                            gsl::span<realtype> dslopesdp) = 0;

    virtual realtype get_value(const realtype t) const = 0;

    virtual realtype get_sensitivity(const realtype t, const int ip) = 0;

    virtual bool get_node_derivative_by_fd() = 0;

    /**
     * @brief Accessor to equidistant_spacing_ member
     * @return equidistant_spacing flag
     */
    bool get_equidistant_spacing() const;

    /**
     * @brief Accessor to logarithmic_parametrization_ member
     * @return logarithmic_parametrization flag
     */
    bool get_logarithmic_parametrization() const;

    int n_nodes() const { return static_cast<int>(nodes_.size()); }

  protected:
    std::vector<realtype> nodes_;

    std::vector<realtype> node_values_;

    std::vector<realtype> coefficients;

    std::vector<realtype> coefficients_extrapolate;

    std::vector<realtype> coefficients_sensi;

    std::vector<realtype> coefficients_extrapolate_sensi;

    virtual void compute_final_value() = 0;

    virtual void compute_final_sensitivity(
      int nplist,
      int spline_offset,
      gsl::span<realtype> dspline_valuesdp,
      gsl::span<realtype> dspline_slopesdp) = 0;

    realtype get_final_value() const;

    void set_final_value(realtype finalValue);

    realtype get_final_sensitivity(const int ip) const;

    void set_final_sensitivity(std::vector<realtype> finalSensitivity);

    /**
     * @brief Switch equidistant spacing of spline nodes on or off
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
    realtype final_value_;

    std::vector<realtype> final_sensitivity_;

    bool equidistant_spacing_ = false;

    bool logarithmic_parametrization_ = false;

}; // class SplineFunction

class HermiteSpline : public AbstractSpline
{
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

    void compute_coefficients() override;

    void compute_coefficients_sensi(int nplist,
                                    int spline_offset,
                                    gsl::span<realtype> dnodesdp,
                                    gsl::span<realtype> dslopesdp) override;

    realtype get_value(const double t) const override;

    realtype get_sensitivity(const double t, const int ip) override;

    /**
     * @brief Accessor to node_derivative_by_FD_ member
     * @return node_derivative_by_FD flag
     */
    bool get_node_derivative_by_fd() override { return node_derivative_by_FD_; }

  private:
    void get_coeffs_sensi_lowlevel(int ip,
                                   int i_node,
                                   int nplist,
                                   int n_spline_coefficients,
                                   int spline_offset,
                                   realtype len,
                                   realtype len_m,
                                   realtype len_p,
                                   gsl::span<realtype> dnodesdp,
                                   gsl::span<realtype> dslopesdp,
                                   gsl::span<realtype> coeffs);

    void handle_inner_derivatives();

    void handle_boundary_conditions();

    void compute_coefficients_extrapolation();

    void compute_coefficients_extrapolation_sensi(
      int nplist,
      int spline_offset,
      gsl::span<realtype> dspline_valuesdp,
      gsl::span<realtype> dspline_slopesdp);

    void compute_final_value() override;

    void compute_final_sensitivity(
      int nplist,
      int spline_offset,
      gsl::span<realtype> dspline_valuesdp,
      gsl::span<realtype> dspline_slopesdp) override;

    std::vector<realtype> node_values_derivative_;

    SplineBoundaryCondition first_bode_bc_ = SplineBoundaryCondition::given;

    SplineBoundaryCondition last_node_bc_ = SplineBoundaryCondition::given;

    SplineExtrapolation first_node_ep_ = SplineExtrapolation::linear;

    SplineExtrapolation last_node_ep_ = SplineExtrapolation::linear;

    bool node_derivative_by_FD_ = false;

}; // class HermiteSpline

} // namespace amici

#endif /* AMICI_SPLINEFUNCTION_H */
