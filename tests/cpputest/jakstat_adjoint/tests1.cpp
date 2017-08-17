#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"
#include <include/amici_hdf5.h>
#include <include/amici_interface_cpp.h>

#include <cstring>
#include "wrapfunctions.h"

TEST_GROUP(groupJakstatAdjoint)
{
    void setup() {

    }

    void teardown() {

    }
};

TEST(groupJakstatAdjoint, testSimulation) {
    Model *model = getModel();

    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/model_jakstat_adjoint/nosensi/options", model);
    ExpData *edata = AMI_HDF5_readSimulationExpData(HDFFILE, udata, "/model_jakstat_adjoint/nosensi/data", model);

    ReturnData *rdata = getSimulationResults(udata, edata);
    CHECK_EQUAL(0, *rdata->status);

    verifyReturnData("/model_jakstat_adjoint/nosensi/results", rdata, udata, model, TEST_ATOL, TEST_RTOL);

    delete model;
    delete edata;
    delete rdata;
    delete udata;
}

TEST(groupJakstatAdjoint, testSensitivityForward) {
    Model *model = getModel();

    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/model_jakstat_adjoint/sensiforward/options", model);
    ExpData *edata = AMI_HDF5_readSimulationExpData(HDFFILE, udata, "/model_jakstat_adjoint/sensiforward/data", model);

    ReturnData *rdata = getSimulationResults(udata, edata);
    CHECK_EQUAL(0, *rdata->status);

    verifyReturnData("/model_jakstat_adjoint/sensiforward/results", rdata, udata, model, TEST_ATOL, TEST_RTOL);

    delete model;
    delete rdata;
    delete edata;
    delete udata;
}

TEST(groupJakstatAdjoint, testSensitivityAdjoint) {
    Model *model = getModel();

    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/model_jakstat_adjoint/sensiadjoint/options", model);
    ExpData *edata = AMI_HDF5_readSimulationExpData(HDFFILE, udata, "/model_jakstat_adjoint/sensiadjoint/data", model);

    ReturnData *rdata = getSimulationResults(udata, edata);
    CHECK_EQUAL(0, *rdata->status);

    verifyReturnData("/model_jakstat_adjoint/sensiadjoint/results", rdata, udata, model, TEST_ATOL, TEST_RTOL);

    delete model;
    delete rdata;
    delete edata;
    delete udata;
}


