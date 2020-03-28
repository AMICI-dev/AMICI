#ifndef AMICI_HDF5_H
#define AMICI_HDF5_H

#include <memory>
#include <string>
#include <vector>

#include <H5Cpp.h>

#include <gsl/gsl-lite.hpp>

#undef AMI_HDF5_H_DEBUG

/* Macros for enabling/disabling  HDF5 error auto-printing
 * AMICI_H5_SAVE_ERROR_HANDLER and AMICI_H5_RESTORE_ERROR_HANDLER must be called
 * within the same context, otherwise the stack handler is lost. */
#define AMICI_H5_SAVE_ERROR_HANDLER                                            \
    herr_t (*old_func)(void *);                                                \
    void *old_client_data;                                                     \
    H5Eget_auto1(&old_func, &old_client_data);                                 \
    H5Eset_auto1(NULL, NULL)

#define AMICI_H5_RESTORE_ERROR_HANDLER H5Eset_auto1(old_func, old_client_data)

namespace amici {

class ReturnData;
class ExpData;
class Model;
class Solver;

namespace hdf5 {

/* Functions for reading and writing AMICI data to/from HDF5 files. */

/**
 * @brief Open the given file for writing. Append if exists, create if not.
 * @param hdf5filename
 * @return
 */
H5::H5File createOrOpenForWriting(std::string const &hdf5filename);

/**
 * @brief Read solver options from HDF5 file
 * @param file hdf5 file handle to read from
 * @param solver solver to set options on
 * @param datasetPath Path inside the HDF5 file
 */
void readSolverSettingsFromHDF5(const H5::H5File &file, Solver &solver,
                                std::string const &datasetPath);

/**
 * @brief Read solver options from HDF5 file
 * @param hdf5Filename Name of HDF5 file
 * @param solver solver to set options on
 * @param datasetPath Path inside the HDF5 file
 */
void readSolverSettingsFromHDF5(std::string const& hdf5Filename, Solver &solver,
                                std::string const &datasetPath);

/**
 * @brief Write solver options from HDF5 file
 * @param hdf5Filenamehdf5 Name of HDF5 file
 * @param solver solver to write options from
 * @param hdf5Location Path inside the HDF5 file
 */
void writeSolverSettingsToHDF5(Solver const& solver,
                              std::string const& hdf5Filename,
                              std::string const& hdf5Location);

/**
 * @brief Write solver options from HDF5 file
 * @param file file handle to read from
 * @param solver solver to write options from
 * @param hdf5Location Path inside the HDF5 file
 */
void writeSolverSettingsToHDF5(Solver const& solver,
                              H5::H5File const& file,
                              std::string const& hdf5Location);

/**
 * @brief Read solver options from HDF5 file
 * @param hdffile Name of HDF5 file
 * @param solver solver to set options on
 * @param datasetPath Path inside the HDF5 file
 */
void readSolverSettingsFromHDF5(std::string const &hdffile, Solver &solver,
                                std::string const &datasetPath);

/**
 * @brief Read model data from HDF5 file
 * @param hdffile Name of HDF5 file
 * @param model model to set data on
 * @param datasetPath Path inside the HDF5 file
 */
void readModelDataFromHDF5(std::string const &hdffile, Model &model,
                           std::string const &datasetPath);

/**
 * @brief Read model data from HDF5 file
 * @param fileId hdf5 file handle to read from
 * @param model model to set data on
 * @param datasetPath Path inside the HDF5 file
 */
void readModelDataFromHDF5(H5::H5File const &file, Model &model,
                           std::string const &datasetPath);

/**
 * @brief Write ReturnData struct to HDF5 dataset
 * @param rdata Data to write
 * @param hdffile Filename of HDF5 file
 * @param datasetPath Full dataset path inside the HDF5 file (will be created)
 */

void writeReturnData(const ReturnData &rdata, H5::H5File const &file,
                     const std::string &hdf5Location);

void writeReturnData(const ReturnData &rdata, std::string const &hdf5Filename,
                     const std::string &hdf5Location);

void writeReturnDataDiagnosis(const ReturnData &rdata, H5::H5File const &file,
                              const std::string &hdf5Location);

/**
 * @brief Create the given group and possibly parents
 * @param file
 * @param groupPath
 * @param recursively
 */
void createGroup(const H5::H5File &file, std::string const &groupPath,
                 bool recursively = true);

/**
 * @brief readSimulationExpData reads AMICI experimental data from attributes in
 * HDF5 file.
 * @param hdf5Filename Name of HDF5 file
 * @param hdf5Root Path inside the HDF5 file to object having ExpData as
 * attributes
 * @param model The model for which data is to be read
 * @return
 */

std::unique_ptr<ExpData> readSimulationExpData(const std::string &hdf5Filename,
                                               const std::string &hdf5Root,
                                               const Model &model);

/**
 * @brief writeSimulationExpData writes AMICI experimental data to attributes in
 * HDF5 file.
 * @param edata The experimental data which is to be written
 * @param hdf5Filename Name of HDF5 file
 * @param hdf5Root Path inside the HDF5 file to object having ExpData as
 * attributes
 */

void writeSimulationExpData(const ExpData &edata, H5::H5File const &file,
                            const std::string &hdf5Location);

/**
 * @brief attributeExists Check whether an attribute with the given name exists
 * on the given dataset
 * @param fileId The HDF5 file object
 * @param datasetPath Dataset of which attributes should be checked
 * @param attributeName Name of the attribute of interest
 * @return
 */
bool attributeExists(H5::H5File const &file, const std::string &optionsObject,
                     const std::string &attributeName);

bool attributeExists(H5::H5Object const &object,
                     const std::string &attributeName);

void createAndWriteInt1DDataset(H5::H5File const &file,
                                std::string const &datasetName,
                                gsl::span<const int> buffer);

void createAndWriteInt2DDataset(H5::H5File const &file,
                                std::string const &datasetName,
                                gsl::span<const int> buffer, hsize_t m,
                                hsize_t n);

void createAndWriteDouble1DDataset(H5::H5File const &file,
                                   std::string const &datasetName,
                                   gsl::span<const double> buffer);

void createAndWriteDouble2DDataset(H5::H5File const &file,
                                   std::string const &datasetName,
                                   gsl::span<const double> buffer, hsize_t m,
                                   hsize_t n);

void createAndWriteDouble3DDataset(H5::H5File const &file,
                                   std::string const &datasetName,
                                   gsl::span<const double> buffer, hsize_t m,
                                   hsize_t n, hsize_t o);

double getDoubleScalarAttribute(const H5::H5File &file,
                                const std::string &optionsObject,
                                const std::string &attributeName);

int getIntScalarAttribute(const H5::H5File &file,
                          const std::string &optionsObject,
                          const std::string &attributeName);

std::vector<int> getIntDataset1D(const H5::H5File &file,
                                 std::string const &name);

std::vector<double> getDoubleDataset1D(const H5::H5File &file,
                                       std::string const &name);

std::vector<double> getDoubleDataset2D(const H5::H5File &file,
                                       std::string const &name, hsize_t &m,
                                       hsize_t &n);

std::vector<double> getDoubleDataset3D(const H5::H5File &file,
                                       std::string const &name, hsize_t &m,
                                       hsize_t &n, hsize_t &o);

/**
 * @brief Check if the given location (group, link or dataset) exists in the
 * given file
 * @param filename
 * @param location
 * @return
 */
bool locationExists(std::string const &filename, std::string const &location);

bool locationExists(H5::H5File const &file, std::string const &location);

} // namespace hdf5
} // namespace amici

#endif
