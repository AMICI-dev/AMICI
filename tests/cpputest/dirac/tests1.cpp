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
    UserData *udata = new UserData();

    udata->qpositivex = new double[model->nx];
    udata->p = new double[model->np];
    udata->k = new double[model->nk];

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
    simulateAndVerifyFromFile(model, "/model_dirac/nosensi/");
    delete model;
}

TEST(groupDirac, testSimulationExpData) {

}

TEST(groupDirac, testSensitivityForward) {
    Model *model = getModel();
    simulateAndVerifyFromFile(model, "/model_dirac/sensiforward/");
    delete model;
}

TEST(groupDirac, testSensitivityState) {

}

TEST(groupDirac, testSensitivityAdjoint) {

}


