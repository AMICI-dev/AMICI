#include "testfunctions.h"
#include <include/amici_hdf5.h>
#include <include/amici_model.h>
#include "include/amici_interface_cpp.h"
#include <cstring>
#include <execinfo.h>
#include <dlfcn.h>    // for dladdr
#include <cxxabi.h>   // for __cxa_demangle
#include <cstdio>
#include <cmath>
#include <unistd.h>

namespace amici {

void simulateAndVerifyFromFile(Model *model, const std::string path)
{
    simulateAndVerifyFromFile(model, HDFFILE, path, TEST_ATOL, TEST_RTOL);
}

void simulateAndVerifyFromFile(Model *model, std::string path, double atol, double rtol)
{
    simulateAndVerifyFromFile(model, HDFFILE, path, atol, rtol);
}


void simulateAndVerifyFromFile(Model *model, const std::string hdffile, std::string path, double atol, double rtol)
{
    // read simulation options
    std::string optionsPath = path + "/options";
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(hdffile.c_str(), optionsPath.c_str(), model);

    std::string measurementPath = path + "/data";
    ExpData *edata = AMI_HDF5_readSimulationExpData(hdffile.c_str(), udata, measurementPath.c_str(), model);

    ReturnData *rdata = getSimulationResults(model, udata, edata);

    std::string resultPath = path + "/results";
    verifyReturnData(hdffile.c_str(), resultPath.c_str(), rdata, udata, model, atol, rtol);

    if(edata)
        delete edata;
    delete rdata;
    delete udata;
}

ExpData *getTestExpData(const UserData *udata, Model *model) {
    ExpData *edata = new ExpData(udata, model);
    return edata;
}

bool withinTolerance(double expected, double actual, double atol, double rtol, int index, const char *name) {
    bool withinTol =  fabs(expected - actual) <= atol || fabs((expected - actual) / (rtol + expected)) <= rtol;

    if(!withinTol && std::isnan(expected) && std::isnan(actual))
        withinTol = true;

    if(!withinTol && std::isinf(expected) && std::isinf(actual))
        withinTol = true;

    if(!withinTol) {
        fprintf(stderr, "ERROR: Expected value %e, but was %e in %s at index %d.\n",expected, actual, name, index);
        fprintf(stderr, "       Relative error: %e (tolerance was %e)\n", fabs((expected - actual) / (rtol + expected)), rtol);
        fprintf(stderr, "       Absolute error: %e (tolerance was %e)\n", fabs(expected - actual), atol);
        printBacktrace(12);
    }

    return withinTol;
}

void checkEqualArray(const double *expected, const double *actual, const int length, double atol, double rtol, const char *name) {
    if(!expected && !actual)
        return;
    if(!length)
        return;

    CHECK_TRUE(expected && actual);

    for(int i = 0; i < length; ++i)
    {
        bool withinTol = withinTolerance(expected[i], actual[i], atol, rtol, i, name);
        CHECK_TRUE(withinTol);
    }
}

void checkEqualArrayStrided(const double *expected, const double *actual, int length, int strideExpected, int strideActual, double atol, double rtol, const char *name) {
    if(!expected && !actual)
        return;

    CHECK_TRUE(expected && actual);

    for(int i = 0; i < length; ++i)
    {
        bool withinTol = withinTolerance(expected[i * strideExpected], actual[i * strideActual], atol, rtol, i, name);
        CHECK_TRUE(withinTol);
    }
}

void verifyReturnData(const char *hdffile, const char* resultPath, const ReturnData *rdata, const UserData *udata, const Model *model, double atol, double rtol) {
    CHECK_FALSE(udata == NULL);
    CHECK_FALSE(rdata == NULL);

    // compare to saved data in hdf file
    hid_t file_id = H5Fopen(hdffile, H5F_ACC_RDONLY, H5P_DEFAULT);

    hsize_t m, n;
    double *expected;

    double statusExp = NAN;
    AMI_HDF5_getDoubleScalarAttribute(file_id, resultPath, "status", &statusExp);
    CHECK_EQUAL((int) statusExp, *rdata->status);
    
    double llhExp = NAN;
    AMI_HDF5_getDoubleScalarAttribute(file_id, resultPath, "llh", &llhExp);
    CHECK_TRUE(withinTolerance(llhExp, *rdata->llh, atol, rtol, 1, "llh"));

    AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "x", &expected, &m, &n);
    checkEqualArray(expected, rdata->x, udata->nt * model->nxtrue, atol, rtol, "x");
    delete[] expected;

//    CHECK_EQUAL(AMICI_O2MODE_FULL, udata->o2mode);

    if(AMI_HDF5_attributeExists(file_id, resultPath, "J")) {
        AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "J", &expected, &m, &n);
        checkEqualArray(expected, rdata->J, model->nx * model->nx, atol, rtol, "J");
        delete[] expected;
    }

    AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "y", &expected, &m, &n);
    checkEqualArray(expected, rdata->y, udata->nt * model->nytrue, atol, rtol, "y");
    delete[] expected;
    
    if(model->nz>0) {
        AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "z", &expected, &m, &n);
        checkEqualArray(expected, rdata->z, udata->nmaxevent * model->nztrue, atol, rtol, "z");
        delete[] expected;
        
        AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "rz", &expected, &m, &n);
        checkEqualArray(expected, rdata->rz, udata->nmaxevent * model->nztrue, atol, rtol, "rz");
        delete[] expected;
        
        AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "sigmaz", &expected, &m, &n);
        checkEqualArray(expected, rdata->sigmaz, udata->nmaxevent * model->nztrue, atol, rtol, "sigmaz");
        delete[] expected;
    }


    AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "xdot", &expected, &m, &n);
    checkEqualArray(expected, rdata->xdot, model->nxtrue, atol, rtol, "xdot");
    delete[] expected;

    if(udata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        verifyReturnDataSensitivities(file_id, resultPath, rdata, udata, model, atol, rtol);
    } else {
        POINTERS_EQUAL(NULL, rdata->sllh);
        POINTERS_EQUAL(NULL, rdata->s2llh);
    }

    H5Fclose(file_id);
}

void verifyReturnDataSensitivities(hid_t file_id, const char* resultPath, const ReturnData *rdata, const UserData*udata, const Model *model, double atol, double rtol) {
    hsize_t m, n, o;
    double *expected;
    int status;

    AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "sllh", &expected, &m, &n);
    checkEqualArray(expected, rdata->sllh, rdata->nplist, atol, rtol, "sllh");
    delete[] expected;

    if(udata->sensi_meth == AMICI_SENSI_FSA) {
        AMI_HDF5_getDoubleArrayAttribute3D(file_id, resultPath, "sx", &expected, &m, &n, &o);
        for(int ip = 0; ip < udata->nplist; ++ip)
            checkEqualArray(&expected[ip * udata->nt * model->nxtrue],
                    &rdata->sx[ip * udata->nt * model->nx],
                    udata->nt * model->nxtrue, atol, rtol, "sx");
        delete[] expected;

        AMI_HDF5_getDoubleArrayAttribute3D(file_id, resultPath, "sy", &expected, &m, &n, &o);
        for(int ip = 0; ip < udata->nplist; ++ip)
            checkEqualArray(&expected[ip * udata->nt * model->nytrue],
                    &rdata->sy[ip * udata->nt * model->ny],
                    udata->nt * model->nytrue, atol, rtol, "sy");
        delete[] expected;

        if(model->nz>0) {
            AMI_HDF5_getDoubleArrayAttribute3D(file_id, resultPath, "sz", &expected, &m, &n, &o);
            for(int ip = 0; ip < udata->nplist; ++ip)
                checkEqualArray(&expected[ip * udata->nmaxevent * model->nztrue],
                                &rdata->sz[ip * udata->nmaxevent * model->nz],
                                udata->nmaxevent * model->nztrue, atol, rtol, "sz");
            delete[] expected;
            
            AMI_HDF5_getDoubleArrayAttribute3D(file_id, resultPath, "srz", &expected, &m, &n, &o);
            for(int ip = 0; ip < udata->nplist; ++ip)
                checkEqualArray(&expected[ip * udata->nmaxevent * model->nztrue],
                                &rdata->srz[ip * udata->nmaxevent * model->nz],
                                udata->nmaxevent * model->nztrue, atol, rtol, "srz");
            delete[] expected;
        }
    }

    status = AMI_HDF5_getDoubleArrayAttribute3D(file_id, resultPath, "ssigmay", &expected, &m, &n, &o);
    if(status == AMICI_SUCCESS) {
        for(int ip = 0; ip < udata->nplist; ++ip)
            checkEqualArray(&expected[ip * udata->nt * model->nytrue],
                    &rdata->ssigmay[ip * udata->nt * model->ny],
                    udata->nt * model->nytrue, atol, rtol, "ssigmay");
        delete[] expected;
    }

    if(model->nz>0) {
        AMI_HDF5_getDoubleArrayAttribute3D(file_id, resultPath, "ssigmaz", &expected, &m, &n, &o);
        for(int ip = 0; ip < udata->nplist; ++ip)
            checkEqualArray(&expected[ip * udata->nmaxevent * model->nztrue],
                            &rdata->ssigmaz[ip * udata->nmaxevent * model->nz],
                            udata->nmaxevent * model->nztrue, atol, rtol, "ssigmaz");
        delete[] expected;
    }

    if(udata->sensi >= AMICI_SENSI_ORDER_SECOND) {
        AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "s2llh", &expected, &m, &n);
        checkEqualArray(expected, rdata->s2llh, (model->nJ-1) * udata->nplist, atol, rtol, "s2llh");
        delete[] expected;
    } else {
        POINTERS_EQUAL(NULL, rdata->s2llh);
        POINTERS_EQUAL(NULL, rdata->s2rz);
    }

}


void printBacktrace(const int nMaxFrames) {
    void *callstack[nMaxFrames];
    char buf[1024];
    int nFrames = backtrace(callstack, nMaxFrames);
    char **symbols = backtrace_symbols(callstack, nFrames);
    
    std::ostringstream trace_buf;
    for (int i = 0; i < nFrames; i++) {
        Dl_info info;
        if (dladdr(callstack[i], &info) && info.dli_sname) {
            char *demangled = NULL;
            int status = -1;
            if (info.dli_sname[0] == '_')
                demangled = abi::__cxa_demangle(info.dli_sname, NULL, 0, &status);
            snprintf(buf, sizeof(buf), "%-3d %*p %s + %zd\n",
                     i, int(2 + sizeof(void*) * 2), callstack[i],
                     status == 0 ? demangled :
                     info.dli_sname == 0 ? symbols[i] : info.dli_sname,
                     (char *)callstack[i] - (char *)info.dli_saddr);
            free(demangled);
        } else {
            snprintf(buf, sizeof(buf), "%-3d %*p %s\n",
                     i, int(2 + sizeof(void*) * 2), callstack[i], symbols[i]);
        }
        trace_buf << buf;
    }
    free(symbols);
    if (nFrames == nMaxFrames)
        trace_buf << "[truncated]\n";
    printf("%s\n",trace_buf.str().c_str());
}

} // namespace amici
