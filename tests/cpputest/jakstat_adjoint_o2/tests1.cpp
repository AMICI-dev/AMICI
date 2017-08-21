#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"
#include <include/amici_hdf5.h>
#include <include/amici_interface_cpp.h>

#include <cstring>
#include "wrapfunctions.h"

TEST_GROUP(groupJakstatAdjointO2)
{
    void setup() {

    }

    void teardown() {

    }
};


TEST(groupJakstatAdjointO2, testSensitivityForward2) {
    Model *model = getModel();

    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/model_jakstat_adjoint/sensi2forward/options", model);
    ExpData *edata = AMI_HDF5_readSimulationExpData(HDFFILE, udata, "/model_jakstat_adjoint/sensi2forward/data", model);
    CHECK_FALSE(edata == NULL);
    CHECK_FALSE(udata == NULL);

    ReturnData *rdata = getSimulationResults(udata, edata);
    CHECK_EQUAL(0, *rdata->status);

    verifyReturnData("/model_jakstat_adjoint/sensi2forward/results", rdata, udata, model, TEST_ATOL, TEST_RTOL);

    delete model;
    delete rdata;
    delete edata;
    delete udata;
}

TEST(groupJakstatAdjointO2, testSensitivityAdjoint2) {
    Model *model = getModel();

    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/model_jakstat_adjoint/sensi2adjoint/options", model);
    ExpData *edata = AMI_HDF5_readSimulationExpData(HDFFILE, udata, "/model_jakstat_adjoint/sensi2adjoint/data", model);
    CHECK_FALSE(edata == NULL);
    CHECK_FALSE(udata == NULL);
    
    ReturnData *rdata = getSimulationResults(udata, edata);
    CHECK_EQUAL(0, *rdata->status);
    
    verifyReturnData("/model_jakstat_adjoint/sensi2adjoint/results", rdata, udata, model, TEST_ATOL, TEST_RTOL);
    
    delete model;
    delete rdata;
    delete edata;
    delete udata;
}




