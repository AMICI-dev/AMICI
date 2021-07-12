#include "testfunctions.h"

#include "wrapfunctions.h"
#include <cstring>

#include <gtest/gtest.h>


TEST(ExampleJakstatAdjointO2, SensitivityForward2)
{
    amici::simulateVerifyWrite("/model_jakstat_adjoint/sensi2forward/");
}

TEST(ExampleJakstatAdjointO2, SensitivityForward2LogParam)
{
    amici::simulateVerifyWrite("/model_jakstat_adjoint/sensi2forwardlogparam/");
}

TEST(ExampleJakstatAdjointO2, SensitivityAdjoint2)
{
    amici::simulateVerifyWrite("/model_jakstat_adjoint/sensi2adjoint/");
}
