#include "testfunctions.h"

#include <cstring>
#include <execinfo.h>
#include <dlfcn.h>    // for dladdr
#include <cxxabi.h>   // for __cxa_demangle
#include <cstdio>
#include <cmath>
#include <unistd.h>
#include <iostream>
#include <format>
#include <span>
#include <string_view>
#include <map>
#include <vector>
#include "gtest/gtest.h"

namespace amici {

namespace generic_model {

extern std::unique_ptr<amici::Model> get_model();

} // namespace generic_model

std::map<std::string, std::vector<std::string_view>> var_names {
    {"p", {"p0", "p1", "p2", "p3", "p4", "p5"}},
    {"k", {"k0", "k1", "k2", "k3", "k4", "p5"}},
    {"x", {"x0", "x1", "x2", "x3", "x4", "x5"}},
    {"y", {"y0", "y1", "y2", "y3", "y4", "y5"}}
};

std::span<std::string_view const> getVariableNames(std::string const& name, int length) {
    auto it = var_names.find(name);
    if(it != var_names.end()) {
        if((int)it->second.size() >= length) {
            return std::span<std::string_view const>(it->second.data(), length);
        } else {
            // need to add more names to var_names above.
            // this is a bit ugly, but we can't generate this dynamically, as
            // we are only passing around references.
            throw std::runtime_error(
                "Variable names for " + name + " with length "
                + std::to_string(length)
                + " not available. Update tests/cpp/testfunctions.cpp:var_names."
            );
        }
    }
    throw std::runtime_error("Variable names for " + name + " with length " + std::to_string(length) + " not found.");
}

void simulateVerifyWrite(const std::string& path)
{
    simulateVerifyWrite(NEW_OPTION_FILE, HDFFILE, HDFFILEWRITE, path, TEST_ATOL, TEST_RTOL);
}

void simulateVerifyWrite(const std::string &path, double atol, double rtol)
{
    simulateVerifyWrite(NEW_OPTION_FILE, HDFFILE, HDFFILEWRITE, path, atol, rtol);
}

void simulateWithDefaultOptions() {
    using namespace amici;
    auto model = amici::generic_model::get_model();
    auto solver = model->create_solver();
    std::unique_ptr<const ExpData> edata;
    auto rdata = run_simulation(*solver, edata.get(), *model);
}

void simulateVerifyWrite(const std::string& hdffileOptions, const std::string& hdffileResults, const std::string& hdffilewrite, const std::string& path, double atol, double rtol)
{
    using namespace amici;
    // read options from file
    std::string optionsPath = path + "/options";
    auto model = amici::generic_model::get_model();
    auto solver = model->create_solver();
    hdf5::read_model_data_from_hdf5(hdffileOptions, *model, optionsPath);
    hdf5::read_solver_settings_from_hdf5(hdffileOptions, *solver, optionsPath);

    // read measurements from file
    std::string measurementPath = path + "/data";
    std::unique_ptr<const ExpData> edata;
    if(hdf5::location_exists(hdffileOptions, measurementPath))
        edata = hdf5::read_exp_data_from_hdf5(hdffileOptions, measurementPath, *model);

    // simulate & verify
    auto rdata = run_simulation(*solver, edata.get(), *model);
    // perform second simulation to check reuse of solver and model object
    auto rdata_resimulation = run_simulation(*solver, edata.get(), *model);
    std::string resultPath = path + "/results";

    // write
    // delete destination group
    H5::H5File in(hdffileOptions.c_str(), H5F_ACC_RDONLY);
    auto out = amici::hdf5::create_or_open_for_writing(hdffilewrite);
    if(hdf5::location_exists(out, path))
        H5Ldelete(out.getId(), path.c_str(), H5P_DEFAULT);
    hdf5::create_group(out, path);
    std::string writePath = path + "/options";
    H5Ocopy(in.getId(), writePath.c_str(), out.getId(), writePath.c_str(), H5P_DEFAULT, H5P_DEFAULT);
    writePath = path + "/data";
    if(edata)
        hdf5::write_exp_data_to_hdf5(*edata, out, writePath);

    writePath = path + "/results";
    hdf5::write_return_data_to_hdf5(*rdata, out, writePath);

    // verify simulated results
    verifyReturnData(hdffileResults, resultPath, rdata.get(), model.get(), atol, rtol);
    // verify resimulated results
    verifyReturnData(hdffileResults, resultPath, rdata_resimulation.get(), model.get(), atol, rtol);
    // verify written results
    verifyReturnData(hdffilewrite, writePath, rdata.get(), model.get(), atol, rtol);
}

std::unique_ptr<ExpData> getTestExpData(Model const& model) {
    return std::unique_ptr<ExpData>(new ExpData(model));
}

bool withinTolerance(double expected, double actual, double atol, double rtol, int index, std::string_view name) {
    bool withinTol =  fabs(expected - actual) <= atol || fabs((expected - actual) / (rtol + expected)) <= rtol;

    if(!withinTol && std::isnan(expected) && std::isnan(actual))
        withinTol = true;

    if(!withinTol && std::isinf(expected) && std::isinf(actual))
        withinTol = true;

    if(!withinTol && (name == "res"  || name == "sres")) {
        // FIXME: Residuals are sign-flipped in old test results.
        // Until the remaining MATLAB tests models are implemented in Python
        // and the correct expected results have been re-generated, we accept
        // both signs.
        withinTol =  fabs(-expected - actual) <= atol || fabs((-expected - actual) / (rtol + -expected)) <= rtol;
    }

    if(!withinTol) {
        fprintf(stderr, "ERROR: Expected value %e, but was %e in %s at index %d.\n",expected, actual, name.data(), index);
        fprintf(stderr, "       Relative error: %e (tolerance was %e)\n", fabs((expected - actual) / (rtol + expected)), rtol);
        fprintf(stderr, "       Absolute error: %e (tolerance was %e)\n", fabs(expected - actual), atol);
        printBacktrace(12);
    }

    return withinTol;
}

void checkEqualArray(std::vector<double> const& expected, std::vector<double> const& actual,
                     double atol, double rtol, std::string const& name) {
    ASSERT_EQ(expected.size(), actual.size());

    for(int i = 0; (unsigned) i < expected.size(); ++i)
    {
        bool withinTol = withinTolerance(expected[i], actual[i], atol, rtol, i, name);
        ASSERT_TRUE(withinTol);
    }
}

void checkEqualArray(const double *expected, const double *actual, const int length, double atol, double rtol, std::string_view name) {
    if(!expected && !actual)
        return;
    if(!length)
        return;

    ASSERT_TRUE(expected && actual);

    for(int i = 0; i < length; ++i)
    {
        bool withinTol = withinTolerance(expected[i], actual[i], atol, rtol, i, name);
        ASSERT_TRUE(withinTol);
    }
}

void checkEqualArrayStrided(const double *expected, const double *actual, int length, int strideExpected, int strideActual, double atol, double rtol, std::string_view name) {
    if(!expected && !actual)
        return;

    ASSERT_TRUE(expected && actual);

    for(int i = 0; i < length; ++i)
    {
        bool withinTol = withinTolerance(expected[i * strideExpected], actual[i * strideActual], atol, rtol, i, name);
        ASSERT_TRUE(withinTol);
    }
}

void verifyReturnData(std::string const& hdffile, std::string const& resultPath,
                      const ReturnData *rdata, const Model *model, double atol, double rtol) {

    ASSERT_FALSE(rdata == nullptr);

    if(!hdf5::location_exists(hdffile, resultPath)) {
        fprintf(stderr, "ERROR: No results available for %s!\n",
                resultPath.c_str());
        return;
    }

    // compare to saved data in hdf file
    H5::H5File file(hdffile.c_str(), H5F_ACC_RDONLY);

    hsize_t m, n;

    std::vector<realtype> expected;

    auto statusExp = hdf5::get_int_scalar_attribute(file, resultPath, "status");
    if(rdata->status != statusExp && !rdata->messages.empty()) {
        std::cerr<<"Messages:"<<std::endl;
        for(const auto& msg: rdata->messages) {
            std::cerr << std::format("[{}][{}] {}", static_cast<int>(msg.severity), msg.identifier, msg.message)<< std::endl;
        }
    }
    ASSERT_EQ(statusExp, rdata->status);

    double llhExp = hdf5::get_double_scalar_attribute(file, resultPath, "llh");
    ASSERT_TRUE(withinTolerance(llhExp, rdata->llh, atol, rtol, 1, "llh"));

    double chi2Exp = hdf5::get_double_scalar_attribute(file, resultPath, "chi2");
    ASSERT_TRUE(withinTolerance(chi2Exp, rdata->chi2, atol, rtol, 1, "chi2"));

    if(hdf5::location_exists(file, resultPath + "/x")) {
        expected = hdf5::get_double_2d_dataset(file, resultPath + "/x", m, n);
        checkEqualArray(expected, rdata->x, atol, rtol, "x");
    } else {
        ASSERT_TRUE(rdata->x.empty());
    }

    //    CHECK_EQUAL(AMICI_O2MODE_FULL, udata->o2mode);

    if(hdf5::location_exists(file, resultPath + "/y")) {
        expected = hdf5::get_double_2d_dataset(file, resultPath + "/y", m, n);
        checkEqualArray(expected, rdata->y, atol, rtol, "y");
    } else {
        ASSERT_TRUE(rdata->y.empty());
    }

    if(hdf5::location_exists(file, resultPath + "/res")) {
        expected = hdf5::get_double_1d_dataset(file, resultPath + "/res");
        checkEqualArray(expected, rdata->res, atol, rtol, "res");
    } else {
        ASSERT_TRUE(rdata->res.empty());
    }

    if(hdf5::location_exists(file, resultPath + "/z")) {
        expected = hdf5::get_double_2d_dataset(file, resultPath + "/z", m, n);
        checkEqualArray(expected, rdata->z, atol, rtol, "z");
    } else {
        ASSERT_TRUE(rdata->z.empty());
    }

    if(hdf5::location_exists(file, resultPath + "/rz")) {
        expected = hdf5::get_double_2d_dataset(file, resultPath + "/rz", m, n);
        checkEqualArray(expected, rdata->rz, atol, rtol, "rz");
    } else {
        ASSERT_TRUE(rdata->rz.empty());
    }

    if(hdf5::location_exists(file, resultPath + "/sigmaz")) {
        expected = hdf5::get_double_2d_dataset(file, resultPath + "/sigmaz", m, n);
        checkEqualArray(expected, rdata->sigmaz, atol, rtol, "sigmaz");
    } else {
        ASSERT_TRUE(rdata->sigmaz.empty());
    }

    expected = hdf5::get_double_1d_dataset(file, resultPath + "/x0");
    checkEqualArray(expected, rdata->x0, atol, rtol, "x0");

    if(rdata->status == AMICI_SUCCESS) {
        // for the failed cases, the stored results don't match
        // since https://github.com/AMICI-dev/AMICI/pull/2349
    expected = hdf5::get_double_1d_dataset(file, resultPath + "/diagnosis/xdot");
    checkEqualArray(expected, rdata->xdot, atol, rtol, "xdot");

    if(hdf5::location_exists(file, resultPath + "/diagnosis/J")) {
        expected = hdf5::get_double_2d_dataset(file, resultPath + "/diagnosis/J", m, n);
        checkEqualArray(expected, rdata->J, atol, rtol, "J");
    } else {
        ASSERT_TRUE(rdata->J.empty());
    }
    }
    if(rdata->sensi >= SensitivityOrder::first) {
        verifyReturnDataSensitivities(file, resultPath, rdata, model, atol, rtol);
    } else {
        ASSERT_EQ(0, rdata->sllh.size());
        ASSERT_EQ(0, rdata->s2llh.size());
    }
}

void verifyReturnDataSensitivities(H5::H5File const& file, std::string const& resultPath,
                                   const ReturnData *rdata, const Model *model, double atol, double rtol) {
    hsize_t m, n, o;
    std::vector<double> expected;
    if(hdf5::location_exists(file, resultPath + "/sllh")) {
        expected = hdf5::get_double_1d_dataset(file, resultPath + "/sllh");
        checkEqualArray(expected, rdata->sllh, atol, rtol, "sllh");
    } else {
        ASSERT_TRUE(rdata->sllh.empty());
    }

    if(rdata->sensi_meth == SensitivityMethod::forward) {

        if(hdf5::location_exists(file, resultPath + "/sx0")) {
            expected = hdf5::get_double_2d_dataset(file, resultPath + "/sx0", m, n);
            checkEqualArray(expected, rdata->sx0, atol, rtol, "sx0");
        } else {
            ASSERT_TRUE(rdata->sx0.empty());
        }

        if(hdf5::location_exists(file, resultPath + "/sres")) {
            expected = hdf5::get_double_2d_dataset(file, resultPath + "/sres", m, n);
            checkEqualArray(expected, rdata->sres, atol, rtol, "sres");
        } else {
            ASSERT_TRUE(rdata->sres.empty());
        }

        if(hdf5::location_exists(file, resultPath + "/FIM")) {
            expected = hdf5::get_double_2d_dataset(file, resultPath + "/FIM", m, n);
            checkEqualArray(expected, rdata->FIM, atol, rtol, "FIM");
        } else {
            ASSERT_TRUE(rdata->FIM.empty());
        }


        /* TODO REMOVE ASAP */
        if(rdata->sensi < SensitivityOrder::second) {
        /* /TODO REMOVE ASAP */

            if(hdf5::location_exists(file, resultPath + "/sx")) {
                expected = hdf5::get_double_3d_dataset(file, resultPath + "/sx", m, n, o);
                for(int ip = 0; ip < model->nplist(); ++ip)
                    checkEqualArray(&expected[ip * model->nt() * model->nxtrue_rdata],
                            &rdata->sx[ip * model->nt() * rdata->nx_rdata],
                            model->nt() * model->nxtrue_rdata, atol, rtol, "sx");
            } else {
                ASSERT_TRUE(rdata->sx.empty());
            }

            if(hdf5::location_exists(file, resultPath + "/sy")) {
                expected = hdf5::get_double_3d_dataset(file, resultPath + "/sy", m, n, o);
                for(int ip = 0; ip < model->nplist(); ++ip)
                    checkEqualArray(&expected[ip * model->nt() * model->nytrue],
                            &rdata->sy[ip * model->nt() * model->ny],
                            model->nt() * model->nytrue, atol, rtol, "sy");
            } else {
                ASSERT_TRUE(rdata->sy.empty());
            }


            if(model->nz>0) {
                expected = hdf5::get_double_3d_dataset(file, resultPath + "/sz", m, n, o);
                for(int ip = 0; ip < model->nplist(); ++ip)
                    checkEqualArray(&expected[ip * model->n_max_event() * model->nztrue],
                            &rdata->sz[ip * model->n_max_event() * model->nz],
                            model->n_max_event() * model->nztrue, atol, rtol, "sz");

                expected = hdf5::get_double_3d_dataset(file, resultPath + "/srz", m, n, o);
                for(int ip = 0; ip < model->nplist(); ++ip)
                    checkEqualArray(&expected[ip * model->n_max_event() * model->nztrue],
                            &rdata->srz[ip * model->n_max_event() * model->nz],
                            model->n_max_event() * model->nztrue, atol, rtol, "srz");
            }

            if(hdf5::location_exists(file, resultPath + "/ssigmay") || !rdata->ssigmay.empty()) {
                expected = hdf5::get_double_3d_dataset(file, resultPath + "/ssigmay", m, n, o);
                for(int ip = 0; ip < model->nplist(); ++ip)
                    checkEqualArray(&expected[ip * model->nt() * model->nytrue],
                            &rdata->ssigmay[ip * model->nt() * model->ny],
                            model->nt() * model->nytrue, atol, rtol, "ssigmay");
            }
            if(hdf5::location_exists(file, resultPath + "/ssigmaz") || !rdata->ssigmaz.empty()) {
                if(model->nz>0) {
                    expected = hdf5::get_double_3d_dataset(file, resultPath + "/ssigmaz", m, n, o);
                    for(int ip = 0; ip < model->nplist(); ++ip)
                        checkEqualArray(&expected[ip * model->n_max_event() * model->nztrue],
                                &rdata->ssigmaz[ip * model->n_max_event() * model->nz],
                                model->n_max_event() * model->nztrue, atol, rtol, "ssigmaz");
                }
            }
        /* TODO REMOVE ASAP */
        }
        /* /TODO REMOVE ASAP */

    }

    if(rdata->sensi >= SensitivityOrder::second) {
        expected = hdf5::get_double_2d_dataset(file, resultPath + "/s2llh", m, n);
        checkEqualArray(expected, rdata->s2llh, atol, rtol, "s2llh");
    } else {
        ASSERT_EQ(0, rdata->s2llh.size());
        ASSERT_EQ(0, rdata->s2rz.size());
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
            char *demangled = nullptr;
            int status = -1;
            if (info.dli_sname[0] == '_')
                demangled = abi::__cxa_demangle(info.dli_sname, nullptr, nullptr, &status);
            snprintf(buf, sizeof(buf), "%-3d %*p %s + %zd\n",
                     i, int(2 + sizeof(void*) * 2), callstack[i],
                     status == 0 ? demangled :
                                   info.dli_sname == nullptr ? symbols[i] : info.dli_sname,
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
