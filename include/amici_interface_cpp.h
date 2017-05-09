#ifndef AMICI_INTERFACE_CPP_H
#define AMICI_INTERFACE_CPP_H

#include "include/amici.h"

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

EXTERNC void initUserDataFields(UserData user_data, ReturnData *rdata, double *pstatus);
EXTERNC  ReturnData *getSimulationResults(UserData *udata, ExpData *edata, int *pstatus);
EXTERNC ReturnData *initReturnData(UserData *udata, int *pstatus);

#endif
