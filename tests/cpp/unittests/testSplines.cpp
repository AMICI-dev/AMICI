#include <amici/splinefunctions.h>

#include <gtest/gtest.h>

#include <cmath>
#include <tuple>
#include <vector>

using std::exp;
using amici::HermiteSpline;
using amici::SplineBoundaryCondition;
using amici::SplineExtrapolation;

#define ASSERT_APPROX(x, x0, rtol) ASSERT_LT(std::abs(((x) - (x0))/(x0)), rtol)

void test_spline_values(HermiteSpline &spline, std::vector<std::tuple<double, double>> &expectations)
{
    for (auto expected : expectations) {
        double time;
        double expected_value;
        std::tie(time, expected_value) = expected;
        ASSERT_DOUBLE_EQ(spline.get_value(time), expected_value);
    }
}

void test_spline_values(HermiteSpline &spline, std::vector<std::tuple<double, double>> &expectations, const double rtol)
{
    for (auto expected : expectations) {
        double time;
        double expected_value;
        std::tie(time, expected_value) = expected;
        ASSERT_APPROX(spline.get_value(time), expected_value, rtol);
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
                         SplineExtrapolation::constant,
                         SplineExtrapolation::constant,
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
}

TEST(Splines, SplineNonUniform)
{
    // Non-uniform grid
    HermiteSpline spline({ 0.0, 0.1, 0.5, 1.0 },
                         { 0.0, 2.0, 0.5, 1.0 },
                         {},
                         SplineBoundaryCondition::given,
                         SplineBoundaryCondition::given,
                         SplineExtrapolation::constant,
                         SplineExtrapolation::constant,
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
}

TEST(Splines, SplineExplicit)
{
    // Derivatives are given explicitly
    HermiteSpline spline({ 0.0, 1.0 },
                         { 0.0, 2.0, 0.5,  1.0, 0.75 },
                         { 1.0, 0.0, 0.1, -0.1, 0.0  },
                         SplineBoundaryCondition::given,
                         SplineBoundaryCondition::given,
                         SplineExtrapolation::constant,
                         SplineExtrapolation::constant,
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
}

TEST(Splines, SplineLogarithmic)
{
    // Logarithmic parametrization
    HermiteSpline spline({ 0.0, 1.0 },
                         { 0.2, 2.0, 0.5, 1.0, 0.75 },
                         {},
                         SplineBoundaryCondition::given,
                         SplineBoundaryCondition::given,
                         SplineExtrapolation::constant,
                         SplineExtrapolation::constant,
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
    test_spline_values(spline, expectations, 1e-14);
}
