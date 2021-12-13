#include <amici/splinefunctions.h>

#include <gtest/gtest.h>

#include <tuple>
#include <vector>

using amici::HermiteSpline;
using amici::SplineBoundaryCondition;
using amici::SplineExtrapolation;

TEST(Splines, test_spline)
{
    HermiteSpline spline({ 1, 2, 3 },
                         { 3, 4, 5 },
                         {},
                         SplineBoundaryCondition::given,
                         SplineBoundaryCondition::given,
                         SplineExtrapolation::constant,
                         SplineExtrapolation::constant,
                         true,
                         false,
                         false);

    // time points and expected values and sensitivities
    std::vector<std::tuple<double, double, int, double>> expectations = {
        // t, expected value, sensitivity parameter idx, expected sensitivity

        // TODO @lcontento
        { 0.0, 1.0, 0, 1.0 }

    };
    for (auto expected : expectations) {
        double time;
        double expected_value;
        int parameter_idx;
        double expected_sensitivity;
        std::tie(time, expected_value, parameter_idx, expected_sensitivity)
            = expected;
        ASSERT_DOUBLE_EQ(spline.get_value(time), expected_value);
        ASSERT_DOUBLE_EQ(spline.get_sensitivity(time, parameter_idx),
                         expected_sensitivity);
    }
}
