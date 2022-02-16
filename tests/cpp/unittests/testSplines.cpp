#include <amici/exception.h>
#include <amici/splinefunctions.h>

#include <gtest/gtest.h>

#include <cmath>
#include <tuple>
#include <vector>

using std::exp;
using amici::HermiteSpline;
using amici::SplineBoundaryCondition;
using amici::SplineExtrapolation;
using amici::AmiException;

#define ASSERT_APPROX(x, x0, rtol, atol) ASSERT_LE(std::abs((x) - (x0)), (atol) + (rtol) * std::abs(x0))

void test_spline_values(HermiteSpline &spline, std::vector<std::tuple<double, double>> &expectations)
{
  for (auto expected : expectations) {
    double time;
    double expected_value;
    std::tie(time, expected_value) = expected;
    ASSERT_DOUBLE_EQ(spline.get_value(time), expected_value);
  }
}

void test_spline_values(HermiteSpline &spline, std::vector<std::tuple<double, double>> &expectations, const double rtol, const double atol)
{
  for (auto expected : expectations) {
    double time;
    double expected_value;
    std::tie(time, expected_value) = expected;
    ASSERT_APPROX(spline.get_value(time), expected_value, rtol, atol);
  }
}

void test_spline_sensitivities(HermiteSpline &spline, std::vector<std::tuple<double, std::vector<double>>> &expectations)
{
  for (auto expected : expectations) {
    double time;
    std::vector<double> expected_values;
    std::tie(time, expected_values) = expected;
    for (std::vector<double>::size_type ip = 0; ip < expected_values.size(); ip++)
      ASSERT_DOUBLE_EQ(spline.get_sensitivity(time, ip), expected_values[ip]);
  }
}

void test_spline_sensitivities(HermiteSpline &spline, std::vector<std::tuple<double, std::vector<double>>> &expectations, const double rtol, const double atol)
{
  for (auto expected : expectations) {
    double time;
    std::vector<double> expected_values;
    std::tie(time, expected_values) = expected;
    for (std::vector<double>::size_type ip = 0; ip < expected_values.size(); ip++)
      ASSERT_APPROX(spline.get_sensitivity(time, ip), expected_values[ip], rtol, atol);
  }
}

TEST(Splines, SplineUniform)
{
    // Uniform grid
    HermiteSpline spline({ 0.0, 1.0 },
                         { 0.0, 2.0, 0.5, 1.0 },
                         {},
                         SplineBoundaryCondition::given,
                         SplineBoundaryCondition::given,
                         SplineExtrapolation::noExtrapolation,
                         SplineExtrapolation::noExtrapolation,
                         true,   // node_derivative_by_FD
                         true,  // equidistant_spacing
                         false); // logarithmic_parametrization

    spline.compute_coefficients();
    std::vector<std::tuple<double, double>> expectations = {
        // t, expected value
        {0.00,  0.0},
        {0.25,  1.74609375},
        {1.0/3, 2.0},
        {0.50,  1.3437499999999996},
        {2.0/3, 0.5},
        {0.75,  0.484375},
        {1.00,  1.0},
    };
    test_spline_values(spline, expectations);
    ASSERT_THROW(spline.get_value(-0.05), AmiException);
    ASSERT_THROW(spline.get_value(1.05), AmiException);
}

TEST(Splines, SplineNonUniform)
{
    // Non-uniform grid
    HermiteSpline spline({ 0.0, 0.1, 0.5, 1.0 },
                         { 0.0, 2.0, 0.5, 1.0 },
                         {},
                         SplineBoundaryCondition::given,
                         SplineBoundaryCondition::given,
                         SplineExtrapolation::noExtrapolation,
                         SplineExtrapolation::noExtrapolation,
                         true,   // node_derivative_by_FD
                         false,  // equidistant_spacing
                         false); // logarithmic_parametrization

    spline.compute_coefficients();
    std::vector<std::tuple<double, double>> expectations = {
        // t, expected value
        {0.00, 0.0},
        {0.05, 1.1484375},
        {0.10, 2.0},
        {0.25, 2.0498046875},
        {0.50, 0.5},
        {0.75, 0.6015625},
        {1.00, 1.0},
    };
    test_spline_values(spline, expectations);
    ASSERT_THROW(spline.get_value(-0.05), AmiException);
    ASSERT_THROW(spline.get_value(1.05), AmiException);
}

TEST(Splines, SplineExplicit)
{
    // Derivatives are given explicitly
    HermiteSpline spline({ 0.0, 1.0 },
                         { 0.0, 2.0, 0.5,  1.0, 0.75 },
                         { 1.0, 0.0, 0.1, -0.1, 0.0  },
                         SplineBoundaryCondition::given,
                         SplineBoundaryCondition::given,
                         SplineExtrapolation::noExtrapolation,
                         SplineExtrapolation::noExtrapolation,
                         false,  // node_derivative_by_FD
                         true,   // equidistant_spacing
                         false); // logarithmic_parametrization

    spline.compute_coefficients();
    std::vector<std::tuple<double, double>> expectations = {
        // t, expected value
        {0.0, 0.0},
        {0.2, 1.8000000000000003},
        {0.25, 2.0},
        {0.4, 1.0243999999999998},
        {0.5, 0.5},
        {0.6, 0.6819999999999999},
        {0.75, 1.0},
        {0.8, 0.9707999999999999},
        {1.0, 0.75},
    };
    test_spline_values(spline, expectations);
    ASSERT_THROW(spline.get_value(-0.05), AmiException);
    ASSERT_THROW(spline.get_value(1.05), AmiException);
}

TEST(Splines, SplineLogarithmic)
{
    // Logarithmic parametrization
    HermiteSpline spline({ 0.0, 1.0 },
                         { 0.2, 2.0, 0.5, 1.0, 0.75 },
                         {},
                         SplineBoundaryCondition::given,
                         SplineBoundaryCondition::given,
                         SplineExtrapolation::noExtrapolation,
                         SplineExtrapolation::noExtrapolation,
                         true,  // node_derivative_by_FD
                         true,  // equidistant_spacing
                         true); // logarithmic_parametrization
    // log-space values [-1.60943791, 0.69314718, -0.69314718, 0, -0.28768207]
    // log-space derivatives [36, 0.3, -4, 0.5, -1.33333333]
    spline.compute_coefficients();
    std::vector<std::tuple<double, double>> expectations = {
        // t, expected value
        {0.0,  0.2},
        {0.2,  2.07939779651678},
        {0.25, 2.0},
        {0.4,  0.947459046694449},
        {0.5,  0.5},
        {0.6,  0.545987404053269},
        {0.75, 1.0},
        {0.8,  0.996753014029391},
        {1.0,  0.75},
    };
    test_spline_values(spline, expectations, 1e-14, 0.0);
    ASSERT_THROW(spline.get_value(-0.05), AmiException);
    ASSERT_THROW(spline.get_value(1.05), AmiException);
}

TEST(Splines, SplineUniformConstantExtrapolation)
{
    // Uniform grid
    HermiteSpline spline({ 0.0, 1.0 },
                         { 0.0, 2.0, 0.5, 1.0 },
                         {},
                         SplineBoundaryCondition::given,
                         SplineBoundaryCondition::given,
                         SplineExtrapolation::constant,
                         SplineExtrapolation::constant,
                         true,   // node_derivative_by_FD
                         true,  // equidistant_spacing
                         false); // logarithmic_parametrization

    spline.compute_coefficients();
    std::vector<std::tuple<double, double>> expectations = {
        // t, expected value
        {-2.00,  0.0},
        {-1.00,  0.0},
        { 0.00,  0.0},
        { 0.25,  1.74609375},
        { 1.0/3, 2.0},
        { 0.50,  1.3437499999999996},
        { 2.0/3, 0.5},
        { 0.75,  0.484375},
        { 1.00,  1.0},
        { 2.00,  1.0},
        { 3.00,  1.0},
    };
    test_spline_values(spline, expectations);
}

TEST(Splines, SplineUniformLinearExtrapolation)
{
    // Uniform grid
    HermiteSpline spline({ 0.0, 1.0 },
                         { 0.0, 2.0, 0.5, 1.0 },
                         {},
                         SplineBoundaryCondition::given,
                         SplineBoundaryCondition::given,
                         SplineExtrapolation::linear,
                         SplineExtrapolation::linear,
                         true,   // node_derivative_by_FD
                         true,  // equidistant_spacing
                         false); // logarithmic_parametrization

    spline.compute_coefficients();
    std::vector<std::tuple<double, double>> expectations = {
        // t, expected value
        {-2.00, -12.0},
        {-1.00,  -6.0},
        { 0.00,   0.0},
        { 0.25,   1.74609375},
        { 1.0/3,  2.0},
        { 0.50,   1.3437499999999996},
        { 2.0/3,  0.5},
        { 0.75,   0.484375},
        { 1.00,   1.0},
        { 2.00,   2.5},
        { 3.00,   4.0},
    };
    test_spline_values(spline, expectations);
}

TEST(Splines, SplineUniformPolynomialExtrapolation)
{
    // Uniform grid
    HermiteSpline spline({ 0.0, 1.0 },
                         { 0.0, 2.0, 0.5, 1.0 },
                         {},
                         SplineBoundaryCondition::given,
                         SplineBoundaryCondition::given,
                         SplineExtrapolation::polynomial,
                         SplineExtrapolation::polynomial,
                         true,   // node_derivative_by_FD
                         true,  // equidistant_spacing
                         false); // logarithmic_parametrization

    spline.compute_coefficients();
    std::vector<std::tuple<double, double>> expectations = {
        // t, expected value
        {-2.00,  429.0},
        {-1.00,   57.0},
        { 0.00,    0.0},
        { 0.25,    1.74609375},
        { 1.0/3,   2.0},
        { 0.50,    1.3437499999999996},
        { 2.0/3,   0.5},
        { 0.75,    0.484375},
        { 1.00,    1.0},
        { 2.00,  -33.5},
        { 3.00, -248.0},
    };
    test_spline_values(spline, expectations);
}

TEST(Splines, SplineUniformPeriodicExtrapolation)
{
    // Uniform grid
    HermiteSpline spline({ 0.0, 1.0 },
                         { 1.0, 2.0, 0.5, 1.0 },
                         {},
                         SplineBoundaryCondition::periodic,
                         SplineBoundaryCondition::periodic,
                         SplineExtrapolation::periodic,
                         SplineExtrapolation::periodic,
                         true,   // node_derivative_by_FD
                         true,  // equidistant_spacing
                         false); // logarithmic_parametrization

    spline.compute_coefficients();
    std::vector<std::tuple<double, double>> expectations = {
        // t, expected value
        {-4.0/3,   0.5},
        {-0.50,    1.2812499999999996},
        { 0.00,    1.0},
        { 0.25,    1.9140625},
        { 1.0/3,   2.0},
        { 0.50,    1.2812499999999996},
        { 2.0/3,   0.5},
        { 0.75,    0.47265625},
        { 1.00,    1.0},
        { 1.25,    1.9140625},
        { 2.75,    0.47265625},
    };
    test_spline_values(spline, expectations);
}

TEST(Splines, SplineNonUniformPeriodicExtrapolation)
{
    // Non-uniform grid
    HermiteSpline spline({ 0.0, 0.1, 0.5, 1.0 },
                         { 1.0, 2.0, 0.5, 1.0 },
                         {},
                         SplineBoundaryCondition::periodic,
                         SplineBoundaryCondition::periodic,
                         SplineExtrapolation::periodic,
                         SplineExtrapolation::periodic,
                         true,   // node_derivative_by_FD
                         false,  // equidistant_spacing
                         false); // logarithmic_parametrization

    spline.compute_coefficients();
    std::vector<std::tuple<double, double>> expectations = {
        // t, expected value
        {-1.90, 2.0},
        {-0.25, 0.3203125},
        { 0.00, 1.0},
        { 0.05, 1.5296875},
        { 0.10, 2.0},
        { 0.25, 1.7568359375},
        { 0.50, 0.5},
        { 0.75, 0.3203125},
        { 1.00, 1.0},
        { 1.50, 0.5},
        { 2.05, 1.5296875},
    };
    test_spline_values(spline, expectations, 1e-14, 0.0);
}

TEST(Splines, SplineUniformSensitivity)
{
    // Uniform grid
    HermiteSpline spline({ 0.0, 1.0 },
                         { 2.5, 3.25, 1.0, 4.5 },
                         {},
                         SplineBoundaryCondition::given,
                         SplineBoundaryCondition::given,
                         SplineExtrapolation::noExtrapolation,
                         SplineExtrapolation::noExtrapolation,
                         true,   // node_derivative_by_FD
                         true,  // equidistant_spacing
                         false); // logarithmic_parametrization
    int n_params = 3;
    std::vector<double> dvaluesdp = {
         3.0,  1.0,  0.0,
         0.0,  0.0,  5.0,
         0.0,  0.0,  0.0,
        -6.0,  1.0,  3.0
    };
    auto dslopesdp = std::vector<double>(spline.n_nodes() * n_params);
    spline.compute_coefficients();
    spline.compute_coefficients_sensi(n_params, 0, dvaluesdp, dslopesdp);
    std::vector<std::tuple<double, std::vector<double>>> expectations = {
        // t, expected values of sensitivities
        {0.00,  {3.0,       1.0,      0.0}},
        {0.25,  {0.539062,  0.179688, 4.45312}},
        {1.0/3, {0.0,       0.0,      5.0}},
        {0.50,  {0.1875,   -0.125,    2.625}},
        {2.0/3, {0.0,       0.0,      0.0}},
        {0.75,  {-1.07812,  0.179688, 0.1875}},
        {1.00,  {-6.0,      1.0,      3.0}},
    };
    test_spline_sensitivities(spline, expectations, 1e-6, 0.0);
    ASSERT_THROW(spline.get_sensitivity(-0.05, 0), AmiException);
    ASSERT_THROW(spline.get_sensitivity( 1.05, 1), AmiException);
}

TEST(Splines, SplineNonUniformSensitivity)
{
    HermiteSpline spline({ 0.0, 0.1, 0.5, 1.0 },
                         { 2.5, 3.25, 1.0, 4.5 },
                         {},
                         SplineBoundaryCondition::given,
                         SplineBoundaryCondition::given,
                         SplineExtrapolation::noExtrapolation,
                         SplineExtrapolation::noExtrapolation,
                         true,   // node_derivative_by_FD
                         false,  // equidistant_spacing
                         false); // logarithmic_parametrization
    int n_params = 3;
    std::vector<double> dvaluesdp = {
         3.0,  1.0,  0.0,
         0.0,  0.0,  5.0,
         0.0,  0.0,  0.0,
        -6.0,  1.0,  3.0
    };
    auto dslopesdp = std::vector<double>(spline.n_nodes() * n_params);
    spline.compute_coefficients();
    spline.compute_coefficients_sensi(n_params, 0, dvaluesdp, dslopesdp);
    std::vector<std::tuple<double, std::vector<double>>> expectations = {
        // t, expected values of sensitivities
        {0.00, { 3.0,     1.0,     0.0}},
        {0.05, { 1.3125,  0.4375,  2.89062}},
        {0.10, { 0.0,     0.0,     5.0}},
        {0.30, {-0.45,   -0.3,     3.6}},
        {0.50, { 0.0,     0.0,     0.0}},
        {0.75, {-2.625,   0.4375,  0.921875}},
        {1.00, {-6.0,     1.0,     3.0}},
    };
    test_spline_sensitivities(spline, expectations, 1e-6, 0.0);
    ASSERT_THROW(spline.get_sensitivity(-0.05, 0), AmiException);
    ASSERT_THROW(spline.get_sensitivity( 1.05, 1), AmiException);
}

TEST(Splines, SplineExplicitSensitivity)
{
    HermiteSpline spline({ 0.0, 1.0 },
                         { 2.5, 3.25, 1.0, 4.5 },
                         { 13.625, 7.5, 1.1585290151921035, 1.0 },
                         SplineBoundaryCondition::given,
                         SplineBoundaryCondition::given,
                         SplineExtrapolation::noExtrapolation,
                         SplineExtrapolation::noExtrapolation,
                         false,  // node_derivative_by_FD
                         true,   // equidistant_spacing
                         false); // logarithmic_parametrization
    int n_params = 3;
    std::vector<double> dvaluesdp = {
         3.0,  1.0,  0.0,
         0.0,  0.0,  5.0,
         0.0,  0.0,  0.0,
        -6.0,  1.0,  3.0
    };
    std::vector<double> dslopesdp = {
         0.0,  0.0,      18.75,
         0.0,  1.0,       3.0,
         4.0, -0.540302,  0.0,
         0.0,  0.0,       0.0,
    };
    std::vector<double>(spline.n_nodes() * n_params);
    spline.compute_coefficients();
    spline.compute_coefficients_sensi(n_params, 0, dvaluesdp, dslopesdp);
    std::vector<std::tuple<double, std::vector<double>>> expectations = {
        // t, expected values of sensitivities
        {0.00,  { 3.0,      1.0,       0.0}},
        {0.25,  { 0.46875,  0.109375,  4.37109}},
        {1.0/3, { 0.0,      0.0,       5.0}},
        {0.50,  {-0.166667, 0.0641793, 2.625}},
        {2.0/3, { 0.0,      0.0,       0.0}},
        {0.75,  {-0.75,     0.130923,  0.46875}},
        {1.00,  {-6.0,      1.0,       3.0}},
    };
    test_spline_sensitivities(spline, expectations, 1e-6, 0.0);
    ASSERT_THROW(spline.get_sensitivity(-0.05, 0), AmiException);
    ASSERT_THROW(spline.get_sensitivity( 1.05, 1), AmiException);
}

TEST(Splines, SplineLogarithmicSensitivity)
{
    HermiteSpline spline({ 0.0, 1.0 },
                         { 2.5, 3.25, 1.0, 4.5 },
                         {},
                         SplineBoundaryCondition::given,
                         SplineBoundaryCondition::given,
                         SplineExtrapolation::noExtrapolation,
                         SplineExtrapolation::noExtrapolation,
                         true,  // node_derivative_by_FD
                         true,  // equidistant_spacing
                         true); // logarithmic_parametrization
    int n_params = 3;
    std::vector<double> dvaluesdp = {
         3.0,  1.0,  0.0,
         0.0,  0.0,  5.0,
         0.0,  0.0,  0.0,
        -6.0,  1.0,  3.0
    };
    auto dslopesdp = std::vector<double>(spline.n_nodes() * n_params);
    spline.compute_coefficients();
    spline.compute_coefficients_sensi(n_params, 0, dvaluesdp, dslopesdp);
    std::vector<std::tuple<double, std::vector<double>>> expectations = {
        // t, expected values of sensitivities
        {0.00, { 3.0,       1.0,        0.0}},
        {0.05, { 1.27483,   0.424942,   6.11077}},
        {0.10, { 0.0,       0.0,        5.0}},
        {0.30, {-0.356731, -0.11891,    1.00886}},
        {0.50, { 0.0,       0.0,        0.0}},
        {0.75, { 0.428375, -0.0713958, -0.214188}},
        {1.00, {-6.0,       1.0,        3.0}},
    };
    test_spline_sensitivities(spline, expectations, 1e-6, 0.0);
    ASSERT_THROW(spline.get_sensitivity(-0.05, 0), AmiException);
    ASSERT_THROW(spline.get_sensitivity( 1.05, 1), AmiException);
}
