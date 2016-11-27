#ifndef AMI_HDF5_H
#define AMI_HDF5_H

#include <hdf5.h>
#include <hdf5_hl.h>
#include <include/amici.h>

#undef AMI_HDF5_H_DEBUG

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

EXTERNC void storeSimulation(const char* fileName, ReturnData *rdata);
EXTERNC UserData *readSimulationUserData(const char* fileName);
EXTERNC ExpData *readSimulationExpData(const char* hdffile, UserData *udata);
EXTERNC void writeReturnData(const char* hdffile, ReturnData *rdata, UserData *udata);
EXTERNC void createAndWriteDouble2DAttribute(hid_t dataset, const char *attributeName, const double *buffer, hsize_t m, hsize_t n);
EXTERNC void createAndWriteDouble3DAttribute(hid_t dataset, const char *attributeName, const double *buffer, hsize_t m, hsize_t n, hsize_t o);
EXTERNC double getDoubleScalarAttribute(hid_t file_id, const char* optionsObject, const char* attributeName);
EXTERNC int getIntScalarAttribute(hid_t file_id, const char* optionsObject, const char* attributeName);
EXTERNC void getDoubleArrayAttribute(hid_t file_id, const char* optionsObject, const char* attributeName, double **destination, hsize_t *length);
EXTERNC void getIntArrayAttribute(hid_t file_id, const char* optionsObject, const char* attributeName, int **destination, hsize_t *length);
EXTERNC void getDoubleArrayAttribute2D(hid_t file_id, const char* optionsObject, const char* attributeName, double **destination, hsize_t *m, hsize_t *n);
EXTERNC void setAttributeIntFromDouble(hid_t file_id, const char *obj_name, const char *attr_name, const double *bufferDouble, size_t size );

#endif
