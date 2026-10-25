#include "wrapfunctions.h"
#include <cstring>
#include <testfunctions.h>

#include <gtest/gtest.h>

TEST(ExampleCalvetti, Simulation)
{
    amici::simulateVerifyWrite("/model_calvetti/nosensi/");
}

TEST(ExampleCalvetti, HasAlgebraicStates)
{
    auto model = amici::generic_model::get_model();
    ASSERT_TRUE(model->has_algebraic_states());
}
