#include "amici/hdf5.h"
#include "amici/amici.h"

#include <hdf5_hl.h>

#include <cassert>
#include <cstring>
#include <algorithm>
#ifdef AMI_HDF5_H_DEBUG
#ifndef __APPLE__
#include <cexecinfo>
#else
#include <execinfo.h>
#endif
#endif
#include <unistd.h>
#include <cmath>


namespace amici {
namespace hdf5 {

/**
 * @brief assertMeasurementDimensionsCompatible
 * @param m
 * @param n
 * @param model
 */
void checkMeasurementDimensionsCompatible(hsize_t m, hsize_t n,
                                          Model const& model) {
    bool compatible = true;
    // if this is rank 1, n and m can be swapped
    if (n == 1) {
        compatible &= (n == (unsigned)model.nt() || n == (unsigned)model.nytrue);
        compatible &= (m == (unsigned)model.nytrue || m == (unsigned)model.nt());
        compatible &= (m * n == (unsigned)model.nytrue * model.nt());
    } else {
        compatible &= (n == (unsigned)model.nytrue);
        compatible &= (m == (unsigned)model.nt());
    }

    if(!compatible)
        throw(AmiException("HDF5 measurement data does not match model. "
                           "Incompatible dimensions."));
}


/**
 * @brief assertEventDimensionsCompatible
 * @param m
 * @param n
 * @param model
 */
void checkEventDimensionsCompatible(hsize_t m, hsize_t n, Model const& model) {
    bool compatible = true;

    // if this is rank 1, n and m can be swapped
    if (n == 1) {
        compatible &= (n == (unsigned)model.nMaxEvent() ||
                       n == (unsigned)model.nztrue);
        compatible &= (m == (unsigned)model.nztrue ||
                       m == (unsigned)model.nMaxEvent());
        compatible &= (m * n == (unsigned)model.nytrue * model.nMaxEvent());
    } else {
        compatible &= (n == (unsigned)model.nztrue);
        compatible &= (m == (unsigned)model.nMaxEvent());
    }

    if(!compatible)
        throw(AmiException("HDF5 event data does not match model. "
                           "Incompatible dimensions."));
}


void createGroup(H5::H5File const& file,
                 std::string const& groupPath,
                 bool recursively) {

    auto groupCreationPropertyList = H5P_DEFAULT;

    if (recursively) {
        groupCreationPropertyList = H5Pcreate(H5P_LINK_CREATE);
        H5Pset_create_intermediate_group(groupCreationPropertyList, 1);
    }

    hid_t group = H5Gcreate(file.getId(), groupPath.c_str(),
                            groupCreationPropertyList,
                            H5P_DEFAULT, H5P_DEFAULT);
    H5Pclose(groupCreationPropertyList);

    if (group < 0)
        throw(AmiException("Failed to create group in hdf5CreateGroup: %s",
                           groupPath.c_str()));
    H5Gclose(group);
}

std::unique_ptr<ExpData> readSimulationExpData(std::string const& hdf5Filename,
                                               std::string const& hdf5Root,
                                               Model const& model) {
    H5::H5File file(hdf5Filename, H5F_ACC_RDONLY);

    hsize_t m, n;

    auto edata = std::unique_ptr<ExpData>(new ExpData(model));

    if (model.ny * model.nt() > 0) {
        if(locationExists(file,  hdf5Root + "/Y")) {
            auto my = getDoubleDataset2D(file, hdf5Root + "/Y", m, n);
            checkMeasurementDimensionsCompatible(m, n, model);
            edata->setObservedData(my);
        } else {
            throw AmiException("Missing %s/Y in %s", hdf5Root.c_str(),
                               hdf5Filename.c_str());
        }

        if(locationExists(file,  hdf5Root + "/Sigma_Y")) {
            auto sigmay = getDoubleDataset2D(file, hdf5Root + "/Sigma_Y", m, n);
            checkMeasurementDimensionsCompatible(m, n, model);
            edata->setObservedDataStdDev(sigmay);
        } else {
            throw AmiException("Missing %s/Sigma_Y in %s", hdf5Root.c_str(),
                               hdf5Filename.c_str());
        }
    }

    if (model.nz * model.nMaxEvent() > 0) {
        if(locationExists(file,  hdf5Root + "/Z")) {
            auto mz = getDoubleDataset2D(file, hdf5Root + "/Z", m, n);
            checkEventDimensionsCompatible(m, n, model);
            edata->setObservedEvents(mz);
        } else {
            throw AmiException("Missing %s/Z in %s", hdf5Root.c_str(),
                               hdf5Filename.c_str());
        }

        if(locationExists(file,  hdf5Root + "/Sigma_Z")) {
            auto sigmaz = getDoubleDataset2D(file, hdf5Root + "/Sigma_Z", m, n);
            checkEventDimensionsCompatible(m, n, model);
            edata->setObservedEventsStdDev(sigmaz);
        } else {
            throw AmiException("Missing %s/Sigma_Z in %s", hdf5Root.c_str(),
                               hdf5Filename.c_str());
        }
    }

    if(locationExists(file,  hdf5Root + "/condition")) {
        edata->fixedParameters = getDoubleDataset1D(file,
                                                    hdf5Root + "/condition");
    }

    if(locationExists(file,  hdf5Root + "/conditionPreequilibration")) {
        edata->fixedParametersPreequilibration = getDoubleDataset1D(
                    file, hdf5Root + "/conditionPreequilibration");
    }

    if(locationExists(file,  hdf5Root + "/conditionPresimulation")) {
        edata->fixedParametersPresimulation = getDoubleDataset1D(
                    file, hdf5Root + "/conditionPresimulation");
    }

    if(attributeExists(file, hdf5Root, "t_presim")) {
        edata->t_presim = getDoubleScalarAttribute(file, hdf5Root, "t_presim");
    }

    if(locationExists(file,  hdf5Root + "/ts")) {
        edata->setTimepoints(getDoubleDataset1D(file, hdf5Root + "/ts"));
    }

    if(attributeExists(file,  hdf5Root,
                        "/reinitializeFixedParameterInitialStates")) {
        edata->reinitializeFixedParameterInitialStates = static_cast<bool>(
            getIntScalarAttribute(file, hdf5Root,
                                  "/reinitializeFixedParameterInitialStates"));
    }

    return edata;
}

void writeSimulationExpData(const ExpData &edata, H5::H5File const& file,
                            const std::string &hdf5Location)
{

    if(!locationExists(file, hdf5Location))
        createGroup(file, hdf5Location);

    if (edata.nt())
        createAndWriteDouble1DDataset(file, hdf5Location + "/ts",
                                      edata.getTimepoints());

    if (!edata.fixedParameters.empty())
        createAndWriteDouble1DDataset(file, hdf5Location + "/condition",
                                      edata.fixedParameters);

    if (!edata.fixedParametersPreequilibration.empty())
        createAndWriteDouble1DDataset(
                    file, hdf5Location + "/conditionPreequilibration",
                    edata.fixedParametersPreequilibration);

    if (!edata.fixedParametersPresimulation.empty())
        createAndWriteDouble1DDataset(
                    file, hdf5Location + "/conditionPresimulation",
                    edata.fixedParametersPresimulation);

    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(), "t_presim",
                             &edata.t_presim, 1);

    if (!edata.getObservedData().empty())
        createAndWriteDouble2DDataset(
                    file, hdf5Location + "/Y", edata.getObservedData(),
                    edata.nt(), edata.nytrue());
    if (!edata.getObservedDataStdDev().empty())
        createAndWriteDouble2DDataset(
                    file, hdf5Location + "/Sigma_Y",
                    edata.getObservedDataStdDev(), edata.nt(), edata.nytrue());
    if (!edata.getObservedEvents().empty())
        createAndWriteDouble2DDataset(file, hdf5Location + "/Z",
                                      edata.getObservedEvents(),
                                      edata.nmaxevent(), edata.nztrue());
    if (!edata.getObservedEventsStdDev().empty())
        createAndWriteDouble2DDataset(file, hdf5Location + "/Sigma_Z",
                                      edata.getObservedEventsStdDev(),
                                      edata.nmaxevent(), edata.nztrue());

    int int_attr = edata.reinitializeFixedParameterInitialStates;
    H5LTset_attribute_int(file.getId(), hdf5Location.c_str(),
                          "reinitializeFixedParameterInitialStates",
                          &int_attr, 1);
}

void writeReturnData(const ReturnData &rdata, H5::H5File const& file, const std::string &hdf5Location)
{

    if(!locationExists(file, hdf5Location))
        createGroup(file, hdf5Location);

    if (!rdata.ts.empty())
        createAndWriteDouble1DDataset(file, hdf5Location + "/t", rdata.ts);

    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(),
                             "llh", &rdata.llh, 1);
    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(),
                             "chi2", &rdata.chi2, 1);

    H5LTset_attribute_int(file.getId(), hdf5Location.c_str(),
                          "status", &rdata.status, 1);

    if (!rdata.sllh.empty())
        createAndWriteDouble1DDataset(file, hdf5Location + "/sllh", rdata.sllh);

    if (!rdata.res.empty())
        createAndWriteDouble1DDataset(file, hdf5Location + "/res", rdata.res);
    if (!rdata.sres.empty())
        createAndWriteDouble2DDataset(file, hdf5Location + "/sres", rdata.sres,
                                      rdata.nt*rdata.nytrue, rdata.nplist);
    if (!rdata.FIM.empty())
        createAndWriteDouble2DDataset(file, hdf5Location + "/FIM",
                                      rdata.FIM, rdata.nplist, rdata.nplist);

    if (!rdata.x0.empty())
        createAndWriteDouble1DDataset(file, hdf5Location + "/x0", rdata.x0);

    if (!rdata.x.empty())
        createAndWriteDouble2DDataset(file, hdf5Location + "/x", rdata.x,
                                      rdata.nt, rdata.nx);

    if (!rdata.y.empty())
        createAndWriteDouble2DDataset(file, hdf5Location + "/y", rdata.y,
                                      rdata.nt, rdata.ny);

    if (!rdata.z.empty())
        createAndWriteDouble2DDataset(file, hdf5Location + "/z", rdata.z,
                                      rdata.nmaxevent, rdata.nz);

    if (!rdata.rz.empty())
        createAndWriteDouble2DDataset(file, hdf5Location + "/rz", rdata.rz,
                                      rdata.nmaxevent, rdata.nz);

    if (!rdata.sigmay.empty())
        createAndWriteDouble2DDataset(file, hdf5Location + "/sigmay",
                                      rdata.sigmay, rdata.nt, rdata.ny);

    if (!rdata.sigmaz.empty())
        createAndWriteDouble2DDataset(file, hdf5Location + "/sigmaz",
                                      rdata.sigmaz, rdata.nmaxevent, rdata.nz);

    if (!rdata.s2llh.empty())
        createAndWriteDouble2DDataset(file, hdf5Location + "/s2llh",
                                      rdata.s2llh, rdata.nJ - 1, rdata.nplist);

    if (!rdata.sx0.empty())
        createAndWriteDouble2DDataset(file, hdf5Location + "/sx0", rdata.sx0,
                                      rdata.nplist, rdata.nx);

    if (!rdata.sx.empty())
        createAndWriteDouble3DDataset(file, hdf5Location + "/sx", rdata.sx,
                                      rdata.nt, rdata.nplist, rdata.nx);

    if (!rdata.sy.empty())
        createAndWriteDouble3DDataset(file, hdf5Location + "/sy", rdata.sy,
                                      rdata.nt, rdata.nplist, rdata.ny);

    if (!rdata.ssigmay.empty())
        createAndWriteDouble3DDataset(file, hdf5Location + "/ssigmay",
                                      rdata.ssigmay, rdata.nt,
                                      rdata.nplist, rdata.ny);

    if (!rdata.sz.empty())
        createAndWriteDouble3DDataset(file, hdf5Location + "/sz", rdata.sz,
                                      rdata.nmaxevent, rdata.nplist, rdata.nz);

    if (!rdata.srz.empty())
        createAndWriteDouble3DDataset(file, hdf5Location + "/srz", rdata.srz,
                                      rdata.nmaxevent, rdata.nplist, rdata.nz);

    if (!rdata.ssigmaz.empty())
        createAndWriteDouble3DDataset(file, hdf5Location + "/ssigmaz",
                                      rdata.ssigmaz,
                                      rdata.nmaxevent, rdata.nplist, rdata.nz);

    writeReturnDataDiagnosis(rdata, file, hdf5Location + "/diagnosis");
}

void writeReturnDataDiagnosis(const ReturnData &rdata,
                              H5::H5File const& file,
                              const std::string& hdf5Location) {

    if(!locationExists(file, hdf5Location))
        createGroup(file, hdf5Location);

    if (!rdata.xdot.empty())
        createAndWriteDouble1DDataset(file, hdf5Location + "/xdot", rdata.xdot);

    if (!rdata.numsteps.empty())
        createAndWriteInt1DDataset(file, hdf5Location + "/numsteps",
                                   rdata.numsteps);

    if (!rdata.numrhsevals.empty())
        createAndWriteInt1DDataset(file, hdf5Location + "/numrhsevals",
                                   rdata.numrhsevals);

    if (!rdata.numerrtestfails.empty())
        createAndWriteInt1DDataset(file, hdf5Location + "/numerrtestfails",
                                   rdata.numerrtestfails);

    if (!rdata.numnonlinsolvconvfails.empty())
        createAndWriteInt1DDataset(file,
                                   hdf5Location + "/numnonlinsolvconvfails",
                                   rdata.numnonlinsolvconvfails);

    if (!rdata.order.empty())
        createAndWriteInt1DDataset(file, hdf5Location + "/order", rdata.order);

    if (!rdata.numstepsB.empty())
        createAndWriteInt1DDataset(file, hdf5Location + "/numstepsB",
                                   rdata.numstepsB);

    if (!rdata.numrhsevalsB.empty())
        createAndWriteInt1DDataset(file, hdf5Location + "/numrhsevalsB",
                                   rdata.numrhsevalsB);

    if (!rdata.numerrtestfailsB.empty())
        createAndWriteInt1DDataset(file, hdf5Location + "/numerrtestfailsB",
                                   rdata.numerrtestfailsB);

    if (!rdata.numnonlinsolvconvfailsB.empty())
        createAndWriteInt1DDataset(
                    file, hdf5Location + "/numnonlinsolvconvfailsB",
                    rdata.numnonlinsolvconvfailsB);

    if (!rdata.preeq_status.empty()) {
        std::vector<int> preeq_status_int (rdata.preeq_status.size());
        for (int i = 0; (unsigned)i < rdata.preeq_status.size(); i++)
            preeq_status_int[i] = static_cast<int>(rdata.preeq_status[i]);
        createAndWriteInt1DDataset(file, hdf5Location + "/preeq_status",
                                   preeq_status_int);
    }

    if (!rdata.preeq_numsteps.empty())
        createAndWriteInt1DDataset(file, hdf5Location + "/preeq_numsteps",
                                   rdata.preeq_numsteps);

    if (!rdata.preeq_numlinsteps.empty())
        createAndWriteInt2DDataset(file, hdf5Location + "/preeq_numlinsteps",
                                   rdata.preeq_numlinsteps,
                                   rdata.newton_maxsteps, 2);

    H5LTset_attribute_int(file.getId(), hdf5Location.c_str(),
                          "preeq_numstepsB", &rdata.preeq_numstepsB, 1);

    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(),
                             "preeq_cpu_time", &rdata.preeq_cpu_time, 1);

    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(),
                             "preeq_cpu_timeB", &rdata.preeq_cpu_timeB, 1);

    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(), "preeq_t",
                             &rdata.preeq_t, 1);

    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(), "preeq_wrms",
                             &rdata.preeq_wrms, 1);

    if (!rdata.posteq_status.empty()) {
        std::vector<int> posteq_status_int (rdata.posteq_status.size());
        for (int i = 0; (unsigned)i < rdata.posteq_status.size(); i++)
            posteq_status_int[i] = static_cast<int>(rdata.posteq_status[i]);
        createAndWriteInt1DDataset(file, hdf5Location + "/posteq_status",
                                   posteq_status_int);
    }

    if (!rdata.posteq_numsteps.empty())
        createAndWriteInt1DDataset(file, hdf5Location + "/posteq_numsteps",
                                   rdata.posteq_numsteps);

    if (!rdata.posteq_numlinsteps.empty())
        createAndWriteInt2DDataset(file, hdf5Location + "/posteq_numlinsteps",
                                   rdata.posteq_numlinsteps,
                                   rdata.newton_maxsteps, 2);

    H5LTset_attribute_int(file.getId(), hdf5Location.c_str(),
                          "posteq_numstepsB", &rdata.posteq_numstepsB, 1);

    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(),
                             "posteq_cpu_time", &rdata.posteq_cpu_time, 1);

    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(),
                             "posteq_cpu_timeB", &rdata.posteq_cpu_timeB, 1);

    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(), "posteq_t",
                             &rdata.posteq_t, 1);

    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(), "posteq_wrms",
                             &rdata.posteq_wrms, 1);

    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(),
                             "cpu_time", &rdata.cpu_time, 1);

    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(),
                             "cpu_timeB", &rdata.cpu_timeB, 1);

    if (!rdata.J.empty())
        createAndWriteDouble2DDataset(file, hdf5Location + "/J", rdata.J,
                                      rdata.nx, rdata.nx);

}


void writeReturnData(ReturnData const& rdata,
                     std::string const& hdf5Filename,
                     std::string const& hdf5Location) {
    auto file = createOrOpenForWriting(hdf5Filename);

    writeReturnData(rdata, file, hdf5Location);
}

double getDoubleScalarAttribute(H5::H5File const& file,
                                std::string const& optionsObject,
                                std::string const& attributeName) {

    double data = NAN;
    herr_t status = H5LTget_attribute_double(file.getId(), optionsObject.c_str(),
                                             attributeName.c_str(), &data);

#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %e\n", attributeName.c_str(), data);
#endif

    if(status < 0)
        throw AmiException("Attribute %s not found for object %s.",
                           attributeName.c_str(), optionsObject.c_str());

    return data;
}

int getIntScalarAttribute(H5::H5File const& file,
                          std::string const& optionsObject,
                          std::string const& attributeName) {
    int data = 0;
    herr_t status = H5LTget_attribute_int(file.getId(), optionsObject.c_str(),
                                          attributeName.c_str(), &data);

#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %d\n", attributeName.c_str(), data);
#endif

    if(status < 0)
        throw AmiException("Attribute %s not found for object %s.",
                           attributeName.c_str(), optionsObject.c_str());

    return data;
}


void createAndWriteInt1DDataset(H5::H5File const& file,
                                std::string const& datasetName,
                                gsl::span<const int> buffer) {
    hsize_t size = buffer.size();
    H5::DataSpace dataspace(1, &size);
    auto dataset = file.createDataSet(datasetName, H5::PredType::NATIVE_INT,
                                      dataspace);
    dataset.write(buffer.data(), H5::PredType::NATIVE_INT);
}

void createAndWriteDouble1DDataset(const H5::H5File &file,
                                   std::string const& datasetName,
                                   gsl::span<const double> buffer) {
    hsize_t size = buffer.size();
    H5::DataSpace dataspace(1, &size);
    auto dataset = file.createDataSet(datasetName, H5::PredType::NATIVE_DOUBLE,
                                      dataspace);
    dataset.write(buffer.data(), H5::PredType::NATIVE_DOUBLE);
}

void createAndWriteDouble2DDataset(const H5::H5File &file,
                                   std::string const& datasetName,
                                   gsl::span<const double> buffer, hsize_t m,
                                   hsize_t n) {
    const hsize_t adims[] {m, n};
    H5::DataSpace dataspace(2, adims);
    auto dataset = file.createDataSet(datasetName, H5::PredType::NATIVE_DOUBLE,
                                      dataspace);
    dataset.write(buffer.data(), H5::PredType::NATIVE_DOUBLE);
}

void createAndWriteInt2DDataset(H5::H5File const& file,
                                std::string const& datasetName,
                                gsl::span<const int> buffer, hsize_t m,
                                hsize_t n) {
    const hsize_t adims[] {m, n};
    H5::DataSpace dataspace(2, adims);
    auto dataset = file.createDataSet(datasetName, H5::PredType::NATIVE_INT,
                                      dataspace);
    dataset.write(buffer.data(), H5::PredType::NATIVE_INT);
}

void createAndWriteDouble3DDataset(H5::H5File const& file,
                                   std::string const& datasetName,
                                   gsl::span<const double> buffer, hsize_t m,
                                   hsize_t n, hsize_t o) {
    const hsize_t adims[] {m, n, o};
    H5::DataSpace dataspace(3, adims);
    auto dataset = file.createDataSet(datasetName, H5::PredType::NATIVE_DOUBLE,
                                      dataspace);
    dataset.write(buffer.data(), H5::PredType::NATIVE_DOUBLE);
}


bool attributeExists(H5::H5File const& file,
                     const std::string &optionsObject,
                     const std::string &attributeName) {
    AMICI_H5_SAVE_ERROR_HANDLER;
    int result = H5Aexists_by_name(file.getId(), optionsObject.c_str(),
                                   attributeName.c_str(), H5P_DEFAULT);
    AMICI_H5_RESTORE_ERROR_HANDLER;
    return result > 0;
}

bool attributeExists(H5::H5Object const& object,
                     const std::string &attributeName) {
    AMICI_H5_SAVE_ERROR_HANDLER;
    int result = H5Aexists(object.getId(), attributeName.c_str());
    AMICI_H5_RESTORE_ERROR_HANDLER;
    return result > 0;
}

void writeSolverSettingsToHDF5(Solver const& solver,
                              std::string const& hdf5Filename,
                              std::string const& hdf5Location) {
    auto file = createOrOpenForWriting(hdf5Filename);

    writeSolverSettingsToHDF5(solver, file, hdf5Location);
}

void writeSolverSettingsToHDF5(Solver const& solver,
                              H5::H5File const& file,
                              const std::string& hdf5Location) {
    if(!locationExists(file, hdf5Location))
        createGroup(file, hdf5Location);

    double dbuffer;
    int ibuffer;

    dbuffer = solver.getAbsoluteTolerance();
    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(),
                             "atol", &dbuffer, 1);

    dbuffer = solver.getRelativeTolerance();
    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(),
                             "rtol", &dbuffer, 1);

    dbuffer = solver.getAbsoluteToleranceFSA();
    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(),
                             "atol_fsa", &dbuffer, 1);

    dbuffer = solver.getRelativeToleranceFSA();
    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(),
                             "rtol_fsa", &dbuffer, 1);

    dbuffer = solver.getAbsoluteToleranceB();
    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(),
                             "atolB", &dbuffer, 1);

    dbuffer = solver.getRelativeToleranceB();
    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(),
                             "rtolB", &dbuffer, 1);

    dbuffer = solver.getAbsoluteToleranceQuadratures();
    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(),
                             "quad_atol", &dbuffer, 1);

    dbuffer = solver.getRelativeToleranceQuadratures();
    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(),
                             "quad_rtol", &dbuffer, 1);

    dbuffer = solver.getAbsoluteToleranceSteadyState();
    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(),
                             "ss_atol", &dbuffer, 1);

    dbuffer = solver.getRelativeToleranceSteadyState();
    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(),
                             "ss_rtol", &dbuffer, 1);

    dbuffer = solver.getAbsoluteToleranceSteadyStateSensi();
    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(),
                             "ss_atol_sensi", &dbuffer, 1);

    dbuffer = solver.getRelativeToleranceSteadyStateSensi();
    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(),
                             "ss_rtol_sensi", &dbuffer, 1);

    ibuffer = static_cast<int>(solver.getMaxSteps());
    H5LTset_attribute_int(file.getId(), hdf5Location.c_str(),
                          "maxsteps", &ibuffer, 1);

    ibuffer = static_cast<int>(solver.getMaxStepsBackwardProblem());
    H5LTset_attribute_int(file.getId(), hdf5Location.c_str(),
                          "maxstepsB", &ibuffer, 1);

    ibuffer = static_cast<int>(solver.getLinearMultistepMethod());
    H5LTset_attribute_int(file.getId(), hdf5Location.c_str(),
                          "lmm", &ibuffer, 1);

    ibuffer = static_cast<int>(solver.getNonlinearSolverIteration());
    H5LTset_attribute_int(file.getId(), hdf5Location.c_str(),
                          "iter", &ibuffer, 1);

    ibuffer = static_cast<int>(solver.getStabilityLimitFlag());
    H5LTset_attribute_int(file.getId(), hdf5Location.c_str(),
                          "stldet", &ibuffer, 1);

    ibuffer = static_cast<int>(solver.getStateOrdering());
    H5LTset_attribute_int(file.getId(), hdf5Location.c_str(),
                          "ordering", &ibuffer, 1);

    ibuffer = static_cast<int>(solver.getInterpolationType());
    H5LTset_attribute_int(file.getId(), hdf5Location.c_str(),
                          "interpType", &ibuffer, 1);

    ibuffer = static_cast<int>(solver.getSensitivityMethod());
    H5LTset_attribute_int(file.getId(), hdf5Location.c_str(),
                          "sensi_meth", &ibuffer, 1);

    ibuffer = static_cast<int>(solver.getSensitivityMethodPreequilibration());
    H5LTset_attribute_int(file.getId(), hdf5Location.c_str(),
                          "sensi_meth_preeq", &ibuffer, 1);

    ibuffer = static_cast<int>(solver.getSensitivityOrder());
    H5LTset_attribute_int(file.getId(), hdf5Location.c_str(),
                          "sensi", &ibuffer, 1);

    ibuffer = static_cast<int>(solver.getNewtonMaxSteps());
    H5LTset_attribute_int(file.getId(), hdf5Location.c_str(),
                          "newton_maxsteps", &ibuffer, 1);

    ibuffer = static_cast<int>(solver.getPreequilibration());
    H5LTset_attribute_int(file.getId(), hdf5Location.c_str(),
                          "newton_preeq", &ibuffer, 1);

    ibuffer = static_cast<int>(solver.getNewtonDampingFactorMode());
    H5LTset_attribute_int(file.getId(), hdf5Location.c_str(),
                          "newton_damping_factor_mode", &ibuffer, 1);

    dbuffer = solver.getNewtonDampingFactorLowerBound();
    H5LTset_attribute_double(file.getId(), hdf5Location.c_str(),
                             "newton_damping_factor_lower_bound", &dbuffer, 1);

    ibuffer = static_cast<int>(solver.getNewtonMaxLinearSteps());
    H5LTset_attribute_int(file.getId(), hdf5Location.c_str(),
                          "newton_maxlinsteps", &ibuffer, 1);

    ibuffer = static_cast<int>(solver.getLinearSolver());
    H5LTset_attribute_int(file.getId(), hdf5Location.c_str(),
                          "linsol", &ibuffer, 1);

    ibuffer = static_cast<int>(solver.getInternalSensitivityMethod());
    H5LTset_attribute_int(file.getId(), hdf5Location.c_str(),
                          "ism", &ibuffer, 1);

    ibuffer = static_cast<int>(solver.getReturnDataReportingMode());
    H5LTset_attribute_int(file.getId(), hdf5Location.c_str(),
                          "rdrm", &ibuffer, 1);
}

void readSolverSettingsFromHDF5(H5::H5File const& file, Solver &solver,
                                const std::string &datasetPath) {

    if(attributeExists(file, datasetPath, "atol")) {
        solver.setAbsoluteTolerance(
                    getDoubleScalarAttribute(file, datasetPath, "atol"));
    }

    if(attributeExists(file, datasetPath, "rtol")) {
        solver.setRelativeTolerance(
                    getDoubleScalarAttribute(file, datasetPath, "rtol"));
    }

    if(attributeExists(file, datasetPath, "atol_fsa")) {
        solver.setAbsoluteToleranceFSA(
                    getDoubleScalarAttribute(file, datasetPath, "atol_fsa"));
    }

    if(attributeExists(file, datasetPath, "rtol_fsa")) {
        solver.setRelativeToleranceFSA(
                    getDoubleScalarAttribute(file, datasetPath, "rtol_fsa"));
    }


    if(attributeExists(file, datasetPath, "atolB")) {
        solver.setAbsoluteToleranceB(
                    getDoubleScalarAttribute(file, datasetPath, "atolB"));
    }

    if(attributeExists(file, datasetPath, "rtolB")) {
        solver.setRelativeToleranceB(
                    getDoubleScalarAttribute(file, datasetPath, "rtolB"));
    }

    if(attributeExists(file, datasetPath, "quad_atol")) {
        solver.setAbsoluteToleranceQuadratures(
                    getDoubleScalarAttribute(file, datasetPath, "quad_atol"));
    }

    if(attributeExists(file, datasetPath, "quad_rtol")) {
        solver.setRelativeToleranceQuadratures(
                    getDoubleScalarAttribute(file, datasetPath, "quad_rtol"));
    }

    if(attributeExists(file, datasetPath, "ss_atol")) {
        solver.setAbsoluteToleranceSteadyState(
                    getDoubleScalarAttribute(file, datasetPath, "ss_atol"));
    }

    if(attributeExists(file, datasetPath, "ss_rtol")) {
        solver.setRelativeToleranceSteadyState(
                    getDoubleScalarAttribute(file, datasetPath, "ss_rtol"));
    }

    if(attributeExists(file, datasetPath, "ss_atol_sensi")) {
        solver.setAbsoluteToleranceSteadyStateSensi(
                    getDoubleScalarAttribute(file, datasetPath,
                                             "ss_atol_sensi"));
    }

    if(attributeExists(file, datasetPath, "ss_rtol_sensi")) {
        solver.setRelativeToleranceSteadyStateSensi(
                    getDoubleScalarAttribute(file, datasetPath,
                                             "ss_rtol_sensi"));
    }

    if(attributeExists(file, datasetPath, "maxsteps")) {
        solver.setMaxSteps(
                    getIntScalarAttribute(file, datasetPath, "maxsteps"));
    }

    if(attributeExists(file, datasetPath, "maxstepsB")) {
        solver.setMaxStepsBackwardProblem(
                    getIntScalarAttribute(file, datasetPath, "maxstepsB"));
    }

    if(attributeExists(file, datasetPath, "lmm")) {
        solver.setLinearMultistepMethod(
                    static_cast<LinearMultistepMethod>(
                        getIntScalarAttribute(file, datasetPath, "lmm")));
    }

    if(attributeExists(file, datasetPath, "iter")) {
        solver.setNonlinearSolverIteration(
                    static_cast<NonlinearSolverIteration>(
                        getIntScalarAttribute(file, datasetPath, "iter")));
    }

    if(attributeExists(file, datasetPath, "stldet")) {
        solver.setStabilityLimitFlag(
                    getIntScalarAttribute(file, datasetPath, "stldet"));
    }

    if(attributeExists(file, datasetPath, "ordering")) {
        solver.setStateOrdering(
                    getIntScalarAttribute(file, datasetPath, "ordering"));
    }

    if(attributeExists(file, datasetPath, "interpType")) {
        solver.setInterpolationType(
                    static_cast<InterpolationType>(
                        getIntScalarAttribute(file, datasetPath,
                                              "interpType")));
    }

    if(attributeExists(file, datasetPath, "sensi_meth")) {
        solver.setSensitivityMethod(
                    static_cast<SensitivityMethod>(
                        getIntScalarAttribute(file, datasetPath, "sensi_meth")));
    }

    if(attributeExists(file, datasetPath, "sensi_meth_preeq")) {
        solver.setSensitivityMethodPreequilibration(
            static_cast<SensitivityMethod>(
                getIntScalarAttribute(file, datasetPath, "sensi_meth_preeq")));
    }

    if(attributeExists(file, datasetPath, "sensi")) {
        solver.setSensitivityOrder(
                    static_cast<SensitivityOrder>(
                        getIntScalarAttribute(file, datasetPath, "sensi")));
    }

    if(attributeExists(file, datasetPath, "newton_maxsteps")) {
        solver.setNewtonMaxSteps(
                    getIntScalarAttribute(file, datasetPath, "newton_maxsteps"));
    }

    if(attributeExists(file, datasetPath, "newton_preeq")) {
        solver.setPreequilibration(
                    getIntScalarAttribute(file, datasetPath, "newton_preeq"));
    }

    if(attributeExists(file, datasetPath, "newton_damping_factor_mode")) {
        solver.setNewtonDampingFactorMode(
                    static_cast<NewtonDampingFactorMode>(
                        getIntScalarAttribute(file, datasetPath, "newton_damping_factor_mode")));
    }

    if(attributeExists(file, datasetPath, "newton_damping_factor_lower_bound")) {
        solver.setNewtonDampingFactorLowerBound(
                    getDoubleScalarAttribute(file, datasetPath, "newton_damping_factor_lower_bound"));
    }

    if(attributeExists(file, datasetPath, "newton_maxlinsteps")) {
        solver.setNewtonMaxLinearSteps(
                    getIntScalarAttribute(file, datasetPath,
                                          "newton_maxlinsteps"));
    }

    if(attributeExists(file, datasetPath, "linsol")) {
        solver.setLinearSolver(
                    static_cast<LinearSolver>(
                        getIntScalarAttribute(file, datasetPath, "linsol")));
    }

    if(attributeExists(file, datasetPath, "ism")) {
        solver.setInternalSensitivityMethod(
                    static_cast<InternalSensitivityMethod>(
                        getIntScalarAttribute(file, datasetPath, "ism")));
    }

    if(attributeExists(file, datasetPath, "rdrm")) {
        solver.setReturnDataReportingMode(
                    static_cast<RDataReporting>(
                        getIntScalarAttribute(file, datasetPath, "rdrm")));
    }
}

void readSolverSettingsFromHDF5(const std::string &hdffile, Solver &solver,
                                const std::string &datasetPath) {
    H5::H5File file(hdffile, H5F_ACC_RDONLY, H5P_DEFAULT);

    readSolverSettingsFromHDF5(file, solver, datasetPath);
}

void readModelDataFromHDF5(const std::string &hdffile, Model &model,
                           const std::string &datasetPath) {
    H5::H5File file(hdffile, H5F_ACC_RDONLY, H5P_DEFAULT);

    readModelDataFromHDF5(file, model, datasetPath);
}

void readModelDataFromHDF5(const H5::H5File &file, Model &model,
                           const std::string &datasetPath) {
    if(attributeExists(file, datasetPath, "tstart")) {
        model.setT0(getDoubleScalarAttribute(file, datasetPath, "tstart"));
    }

    if(locationExists(file, datasetPath + "/pscale")) {
        auto pscaleInt = getIntDataset1D(file, datasetPath + "/pscale");
        std::vector<ParameterScaling> pscale(pscaleInt.size());
        for(int i = 0; (unsigned)i < pscaleInt.size(); ++i)
            pscale[i] = static_cast<ParameterScaling>(pscaleInt[i]);
        model.setParameterScale(pscale);
    } else if (attributeExists(file, datasetPath, "pscale")) {
        // if pscale is the same for all parameters,
        // it can be set as scalar attribute for convenience
        model.setParameterScale(
                    static_cast<ParameterScaling>(
                        getDoubleScalarAttribute(file,  datasetPath, "pscale")));
    }

    if(attributeExists(file, datasetPath, "nmaxevent")) {
        model.setNMaxEvent(getIntScalarAttribute(file, datasetPath, "nmaxevent"));
    }

    if(attributeExists(file, datasetPath, "steadyStateSensitivityMode")) {
        model.setSteadyStateSensitivityMode(
                    static_cast<SteadyStateSensitivityMode>(
                        getIntScalarAttribute(file, datasetPath,
                                              "steadyStateSensitivityMode")));
    }

    if(locationExists(file, datasetPath + "/theta")) {
        model.setParameters(getDoubleDataset1D(file, datasetPath + "/theta"));
    }

    if(locationExists(file, datasetPath + "/kappa")) {
        model.setFixedParameters(getDoubleDataset1D(file, datasetPath + "/kappa"));
    }

    if(locationExists(file, datasetPath + "/ts")) {
        model.setTimepoints(getDoubleDataset1D(file, datasetPath + "/ts"));
    }

    if(locationExists(file, datasetPath + "/sens_ind")) {
        auto sensInd = getIntDataset1D(file, datasetPath + "/sens_ind");
        model.setParameterList(sensInd);
    }

    if(locationExists(file, datasetPath + "/x0")) {
        auto x0 = getDoubleDataset1D(file, datasetPath + "/x0");
        if(!x0.empty())
            model.setInitialStates(x0);
    }

    if(locationExists(file, datasetPath + "/sx0")) {
        hsize_t length0 = 0;
        hsize_t length1 = 0;
        auto sx0 = getDoubleDataset2D(file, datasetPath + "/sx0",
                                      length0, length1);
        if(!sx0.empty()) {
            if (length0 != (unsigned) model.nplist()
                    && length1 != (unsigned) model.nx_rdata)
                throw(AmiException("Dimension mismatch when reading sx0. "
                                   "Expected %dx%d, got %llu, %llu.",
                                   model.nx_rdata, model.nplist(), length0, length1));
            model.setUnscaledInitialStateSensitivities(sx0);
        }
    }

}

H5::H5File createOrOpenForWriting(const std::string &hdf5filename)
{
    AMICI_H5_SAVE_ERROR_HANDLER;
    try {
        H5::H5File file(hdf5filename, H5F_ACC_RDWR);
        AMICI_H5_RESTORE_ERROR_HANDLER;
        return file;
    } catch(...) {
        AMICI_H5_RESTORE_ERROR_HANDLER;
        return H5::H5File(hdf5filename, H5F_ACC_EXCL);
    }
}

bool locationExists(const H5::H5File &file, const std::string &location)
{
    AMICI_H5_SAVE_ERROR_HANDLER;
    auto result = H5Lexists(file.getId(), location.c_str(), H5P_DEFAULT) > 0;
    AMICI_H5_RESTORE_ERROR_HANDLER;
    return result;
}

bool locationExists(const std::string &filename, const std::string &location)
{
    H5::H5File file(filename, H5F_ACC_RDONLY);
    return locationExists(file, location);
}

std::vector<int> getIntDataset1D(const H5::H5File &file,
                                 std::string const& name) {
    auto dataset = file.openDataSet(name);
    auto dataspace = dataset.getSpace();

    int rank = dataspace.getSimpleExtentNdims();
    if(rank != 1)
        throw(AmiException("Expected array of rank 1 in %s", name.c_str()));

    hsize_t dim;
    dataspace.getSimpleExtentDims(&dim);
    std::vector<int> result(dim);
    if(!result.empty())
        dataset.read(result.data(), H5::PredType::NATIVE_INT);
    return result;
}


std::vector<double> getDoubleDataset1D(const H5::H5File &file,
                                       const std::string &name)
{
    auto dataset = file.openDataSet(name);
    auto dataspace = dataset.getSpace();

    int rank = dataspace.getSimpleExtentNdims();
    if(rank != 1)
        throw(AmiException("Expected array of rank 1 in %s", name.c_str()));

    hsize_t dim;
    dataspace.getSimpleExtentDims(&dim);
    std::vector<double> result(dim);
    if(!result.empty())
        dataset.read(result.data(), H5::PredType::NATIVE_DOUBLE);

    return result;

}

std::vector<double> getDoubleDataset2D(const H5::H5File &file,
                                       const std::string &name,
                                       hsize_t &m, hsize_t &n)
{
    m = n = 0;

    auto dataset = file.openDataSet(name);
    auto dataspace = dataset.getSpace();

    int rank = dataspace.getSimpleExtentNdims();
    if(rank != 2)
        throw(AmiException("Expected array of rank 2 in %s", name.c_str()));

    hsize_t dims[2];
    dataspace.getSimpleExtentDims(dims);
    m = dims[0];
    n = dims[1];

    std::vector<double> result(m * n);
    if(!result.empty())
        dataset.read(result.data(), H5::PredType::NATIVE_DOUBLE);

    return result;
}

std::vector<double> getDoubleDataset3D(const H5::H5File &file,
                                       const std::string &name,
                                       hsize_t &m, hsize_t &n, hsize_t &o)
{
    m = n = o = 0;

    auto dataset = file.openDataSet(name);
    auto dataspace = dataset.getSpace();

    int rank = dataspace.getSimpleExtentNdims();
    if(rank != 3)
        throw(AmiException("Expected array of rank 3 in %s", name.c_str()));

    hsize_t dims[3];
    dataspace.getSimpleExtentDims(dims);
    m = dims[0];
    n = dims[1];
    o = dims[2];

    std::vector<double> result(m * n * o);
    if(!result.empty())
        dataset.read(result.data(), H5::PredType::NATIVE_DOUBLE);

    return result;
}

} // namespace hdf5
} // namespace amici
