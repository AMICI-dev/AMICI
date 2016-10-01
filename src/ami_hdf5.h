#ifndef AMI_HDF5_H
#define AMI_HDF5_H

#include <hdf5.h>
#include <hdf5_hl.h>
#include <include/amici.h>

#undef AMI_HDF5_H_DEBUG

void storeSimulation(const char* fileName, ReturnData *rdata);
UserData *readSimulationUserData(const char* fileName);
ExpData *readSimulationExpData(const char* hdffile, UserData *udata);
void writeReturnData(const char* hdffile, ReturnData *rdata, UserData *udata);
void createAndWriteDouble2DAttribute(hid_t dataset, const char *attributeName, const double *buffer, hsize_t m, hsize_t n);
void createAndWriteDouble3DAttribute(hid_t dataset, const char *attributeName, const double *buffer, hsize_t m, hsize_t n, hsize_t o);
double getDoubleScalarAttribute(hid_t file_id, const char* optionsObject, const char* attributeName);
int getIntScalarAttribute(hid_t file_id, const char* optionsObject, const char* attributeName);
void getDoubleArrayAttribute(hid_t file_id, const char* optionsObject, const char* attributeName, double **destination, hsize_t *length);
void getIntArrayAttribute(hid_t file_id, const char* optionsObject, const char* attributeName, int **destination, hsize_t *length);
void getDoubleArrayAttribute2D(hid_t file_id, const char* optionsObject, const char* attributeName, double **destination, hsize_t *m, hsize_t *n);
void setAttributeIntFromDouble(hid_t file_id, const char *obj_name, const char *attr_name, const double *bufferDouble, size_t size );

#endif
