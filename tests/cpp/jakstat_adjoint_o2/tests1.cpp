#include "testfunctions.h"

#include "wrapfunctions.h"
#include <cstring>

#include <gtest/gtest.h>


TEST(groupJakstatAdjointO2, testSensitivityForward2)
{
    amici::simulateVerifyWrite("/model_jakstat_adjoint/sensi2forward/");
}

TEST(groupJakstatAdjointO2, testSensitivityForward2LogParam)
{
    amici::simulateVerifyWrite("/model_jakstat_adjoint/sensi2forwardlogparam/");
}

TEST(groupJakstatAdjointO2, testSensitivityAdjoint2)
{
    amici::simulateVerifyWrite("/model_jakstat_adjoint/sensi2adjoint/");
}
