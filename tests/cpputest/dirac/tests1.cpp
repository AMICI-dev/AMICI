#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include <testfunctions.h>
#include <include/amici_hdf5.h>
#include <include/amici_interface_cpp.h>
#include <cstring>
#include "wrapfunctions.h"

TEST_GROUP(groupDirac)
{
    void setup() {

    }

    void teardown() {

    }
};


UserData *getTestUserData() {
    Model *model = getModel();
    UserData *udata = new UserData(getUserData());

    udata->qpositivex = new double[model->nx];
    udata->p = new double[model->np];
    udata->k = new double[model->nk];
    udata->idlist = new realtype[model->np]();
    for(int i = 0; i < model->np; ++i)
        udata->idlist[i] = 0;

    delete model;
    return udata;
}


/*
 * Test for mem leaks in UserData initalization / destruction
 */
TEST(groupDirac, testCreateAndFreeUserData) {
    UserData *udata = getTestUserData();

    delete udata;
}

/*
 * Test for mem leaks in ExpData initalization / destruction
 */

TEST(groupDirac, testCreateAndFreeExpData) {
    UserData udata;
    
    Model *model = getModel();

    ExpData *edata = getTestExpData(&udata, model);
    
    delete model;
    delete edata;
}

/*
 * Test for mem leaks in ReturnData initalization / destruction
 */

TEST(groupDirac, testCreateAndFreeReturnData) {
    Model *model = getModel();

    UserData udata;
    ReturnData *rdata = new ReturnData(&udata, model);
    delete model;

    delete rdata;
}


TEST(groupDirac, testSimulation) {
    Model *model = getModel();

    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/model_dirac/nosensi/options", model);
    ExpData *edata = NULL;

    ReturnData *rdata = getSimulationResults(udata, edata);
    CHECK_EQUAL(0, *rdata->status);

    // TODO proper paths /testDirac1/...
    verifyReturnData("/model_dirac/nosensi/results", rdata, udata, model, TEST_ATOL, TEST_RTOL);

    delete model;
    delete rdata;
    delete udata;
}

TEST(groupDirac, testSimulationExpData) {

}

TEST(groupDirac, testSensitivityForward) {
    Model *model = getModel();

    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/model_dirac/sensiforward/options", model);
    ExpData *edata = NULL;

    ReturnData *rdata = getSimulationResults(udata, edata);
    CHECK_EQUAL(0, *rdata->status);

    // TODO proper paths /testDirac1/...
    verifyReturnData("/model_dirac/sensiforward/results", rdata, udata, model, TEST_ATOL, TEST_RTOL);

    delete model;
    delete rdata;
    delete udata;
}

TEST(groupDirac, testSensitivityState) {

}

TEST(groupDirac, testSensitivityAdjoint) {

}


