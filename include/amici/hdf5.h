#ifndef AMICI_HDF5_H
#define AMICI_HDF5_H

#include <memory>
#include <string>
#include <vector>

#include <H5Cpp.h>

#include <gsl/gsl-lite.hpp>

/* Macros for enabling/disabling  HDF5 error auto-printing
 * AMICI_H5_SAVE_ERROR_HANDLER and AMICI_H5_RESTORE_ERROR_HANDLER must be called
 * within the same context, otherwise the stack handler is lost. */
#define AMICI_H5_SAVE_ERROR_HANDLER                                            \
    herr_t (*old_func)(void*);                                                 \
    void* old_client_data;                                                     \
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
 * @brief Open the given file for writing.
 *
 * Append if exists, create if not.
 * @param hdf5filename File to open
 * @return File object
 */
H5::H5File createOrOpenForWriting(std::string const& hdf5filename);

/**
 * @brief Read solver options from HDF5 file.
 * @param file HDF5 file to read from
 * @param solver Solver to set options on
 * @param datasetPath Path inside the HDF5 file
 */
void readSolverSettingsFromHDF5(
    const H5::H5File& file, Solver& solver, std::string const& datasetPath
);

/**
 * @brief Write solver options to HDF5 file.
 * @param hdf5Filename Name of HDF5 file to write to
 * @param solver Solver to write options from
 * @param hdf5Location Path inside the HDF5 file
 */
void writeSolverSettingsToHDF5(
    Solver const& solver, std::string const& hdf5Filename,
    std::string const& hdf5Location
);

/**
 * @brief Write solver options to HDF5 file.
 * @param file File to read from
 * @param solver Solver to write options from
 * @param hdf5Location Path inside the HDF5 file
 */
void writeSolverSettingsToHDF5(
    Solver const& solver, H5::H5File const& file,
    std::string const& hdf5Location
);

/**
 * @brief Read solver options from HDF5 file.
 * @param hdffile Name of HDF5 file
 * @param solver Solver to set options on
 * @param datasetPath Path inside the HDF5 file
 */
void readSolverSettingsFromHDF5(
    std::string const& hdffile, Solver& solver, std::string const& datasetPath
);

/**
 * @brief Read model data from HDF5 file.
 * @param hdffile Name of HDF5 file
 * @param model Model to set data on
 * @param datasetPath Path inside the HDF5 file
 */
void readModelDataFromHDF5(
    std::string const& hdffile, Model& model, std::string const& datasetPath
);

/**
 * @brief Read model data from HDF5 file.
 * @param file HDF5 file handle to read from
 * @param model Model to set data on
 * @param datasetPath Path inside the HDF5 file
 */
void readModelDataFromHDF5(
    H5::H5File const& file, Model& model, std::string const& datasetPath
);

/**
 * @brief Write ReturnData to HDF5 file.
 * @param rdata Data to write
 * @param file HDF5 file to write to
 * @param hdf5Location Full dataset path inside the HDF5 file (will be created)
 */

void writeReturnData(
    ReturnData const& rdata, H5::H5File const& file,
    std::string const& hdf5Location
);

/**
 * @brief Write ReturnData to HDF5 file.
 * @param rdata Data to write
 * @param hdf5Filename Filename of HDF5 file
 * @param hdf5Location Full dataset path inside the HDF5 file (will be created)
 */

void writeReturnData(
    ReturnData const& rdata, std::string const& hdf5Filename,
    std::string const& hdf5Location
);

/**
 * @brief Write ReturnData diagnosis data to HDF5 file.
 * @param rdata Data to write
 * @param file HDF5 file to write to
 * @param hdf5Location Full dataset path inside the HDF5 file (will be created)
 */
void writeReturnDataDiagnosis(
    ReturnData const& rdata, H5::H5File const& file,
    std::string const& hdf5Location
);

/**
 * @brief Create the given group and possibly parents.
 * @param file HDF5 file to write to
 * @param groupPath Path to the group to be created
 * @param recursively Create intermediary groups
 */
void createGroup(
    const H5::H5File& file, std::string const& groupPath,
    bool recursively = true
);

/**
 * @brief Read AMICI ExpData data from HDF5 file.
 * @param hdf5Filename Name of HDF5 file
 * @param hdf5Root Path inside the HDF5 file to object having ExpData
 * @param model The model for which data is to be read
 * @return ExpData created from data in the given location
 */

std::unique_ptr<ExpData> readSimulationExpData(
    std::string const& hdf5Filename, std::string const& hdf5Root,
    Model const& model
);

/**
 * @brief Write AMICI experimental data to HDF5 file.
 * @param edata The experimental data which is to be written
 * @param file Name of HDF5 file
 * @param hdf5Location Path inside the HDF5 file to object having ExpData
 */

void writeSimulationExpData(
    ExpData const& edata, H5::H5File const& file,
    std::string const& hdf5Location
);

/**
 * @brief Check whether an attribute with the given name exists
 * on the given dataset.
 * @param file The HDF5 file object
 * @param optionsObject Dataset of which attributes should be checked
 * @param attributeName Name of the attribute of interest
 * @return `true` if attribute exists, `false` otherwise
 */
bool attributeExists(
    H5::H5File const& file, std::string const& optionsObject,
    std::string const& attributeName
);

/**
 * @brief Check whether an attribute with the given name exists
 * on the given object.
 * @param object An HDF5 object
 * @param attributeName Name of the attribute of interest
 * @return `true` if attribute exists, `false` otherwise
 */
bool attributeExists(
    H5::H5Object const& object, std::string const& attributeName
);

/**
 * @brief Create and write to 1-dimensional native integer dataset.
 * @param file HDF5 file object
 * @param datasetName Name of dataset to create
 * @param buffer Data to write to dataset
 */
void createAndWriteInt1DDataset(
    H5::H5File const& file, std::string const& datasetName,
    gsl::span<int const> buffer
);

/**
 * @brief Create and write to 2-dimensional native integer dataset.
 * @param file HDF5 file object
 * @param datasetName Name of dataset to create
 * @param buffer Flattened data to write to dataset (assuming row-major)
 * @param m Number of rows in buffer
 * @param n Number of columns buffer
 */
void createAndWriteInt2DDataset(
    H5::H5File const& file, std::string const& datasetName,
    gsl::span<int const> buffer, hsize_t m, hsize_t n
);

/**
 * @brief Create and write to 1-dimensional native double dataset.
 * @param file HDF5 file object
 * @param datasetName Name of dataset to create
 * @param buffer Data to write to dataset
 */
void createAndWriteDouble1DDataset(
    H5::H5File const& file, std::string const& datasetName,
    gsl::span<double const> buffer
);

/**
 * @brief Create and write to 2-dimensional native double dataset.
 * @param file HDF5 file object
 * @param datasetName Name of dataset to create
 * @param buffer Flattened data to write to dataset (assuming row-major)
 * @param m Number of rows in buffer
 * @param n Number of columns buffer
 */

void createAndWriteDouble2DDataset(
    H5::H5File const& file, std::string const& datasetName,
    gsl::span<double const> buffer, hsize_t m, hsize_t n
);

/**
 * @brief Create and write to 3-dimensional native double dataset.
 * @param file HDF5 file object
 * @param datasetName Name of dataset to create
 * @param buffer Flattened data to write to dataset (assuming row-major)
 * @param m Length of first dimension in buffer
 * @param n Length of first dimension in buffer
 * @param o Length of first dimension in buffer
 */

void createAndWriteDouble3DDataset(
    H5::H5File const& file, std::string const& datasetName,
    gsl::span<double const> buffer, hsize_t m, hsize_t n, hsize_t o
);

/**
 * @brief Read string attribute from HDF5 object.
 * @param file HDF5 file
 * @param optionsObject Object to read attribute from
 * @param attributeName Name of attribute to read
 * @return Attribute value
 */
std::string getStringAttribute(
    H5::H5File const& file, std::string const& optionsObject,
    std::string const& attributeName
);

/**
 * @brief Read scalar native double attribute from HDF5 object.
 * @param file HDF5 file
 * @param optionsObject Object to read attribute from
 * @param attributeName Name of attribute to read
 * @return Attribute value
 */
double getDoubleScalarAttribute(
    const H5::H5File& file, std::string const& optionsObject,
    std::string const& attributeName
);

/**
 * @brief Read scalar native integer attribute from HDF5 object.
 * @param file HDF5 file
 * @param optionsObject Object to read attribute from
 * @param attributeName Name of attribute to read
 * @return Attribute value
 */

int getIntScalarAttribute(
    const H5::H5File& file, std::string const& optionsObject,
    std::string const& attributeName
);

/**
 * @brief Read 1-dimensional native integer dataset from HDF5 file.
 * @param file HDF5 file object
 * @param name Name of dataset to read
 * @return Data read
 */
std::vector<int>
getIntDataset1D(const H5::H5File& file, std::string const& name);

/**
 * @brief Read 1-dimensional native double dataset from HDF5 file.
 * @param file HDF5 file object
 * @param name Name of dataset to read
 * @return Data read
 */

std::vector<double>
getDoubleDataset1D(const H5::H5File& file, std::string const& name);

/**
 * @brief Read 2-dimensional native double dataset from HDF5 file.
 * @param file HDF5 file object
 * @param name Name of dataset to read
 * @param m Number of rows in the dataset
 * @param n Number of columns in the dataset
 * @return Flattened data (row-major)
 */

std::vector<double> getDoubleDataset2D(
    const H5::H5File& file, std::string const& name, hsize_t& m, hsize_t& n
);

/**
 * @brief Read 3-dimensional native double dataset from HDF5 file.
 * @param file HDF5 file object
 * @param name Name of dataset to read
 * @param m Length of first dimension in dataset
 * @param n Length of first dimension in dataset
 * @param o Length of first dimension in dataset
 * @return Flattened data (row-major)
 */

std::vector<double> getDoubleDataset3D(
    const H5::H5File& file, std::string const& name, hsize_t& m, hsize_t& n,
    hsize_t& o
);

/**
 * @brief Check if the given location (group, link or dataset) exists in the
 * given file.
 * @param filename HDF5 filename
 * @param location Location to test for
 * @return `true` if exists, `false` otherwise
 */
bool locationExists(std::string const& filename, std::string const& location);

/**
 * @brief Check if the given location (group, link or dataset) exists in the
 * given file.
 * @param file HDF5 file object
 * @param location Location to test for
 * @return `true` if exists, `false` otherwise
 */

bool locationExists(H5::H5File const& file, std::string const& location);

} // namespace hdf5
} // namespace amici

#endif
