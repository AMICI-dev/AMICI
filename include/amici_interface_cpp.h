#ifndef AMICI_INTERFACE_CPP_H
#define AMICI_INTERFACE_CPP_H

#include "amici_defines.h"
#include "include/amici_misc.h"
#include <include/edata.h>
#include <include/rdata.h>
#include <include/udata.h>

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

/**
 * getSimulationResults runs the forward an backwards simulation and returns
 * results in a ReturnData struct
 *
 * @param[in] udata pointer to the user data struct @type UserData
 * @param[in] edata pointer to the experimental data struct @type ExpData
 * @return rdata data struct with simulation results @type ReturnData
 */

EXTERNC ReturnData *getSimulationResults(Model *model, UserData *udata,
                                         const ExpData *edata);

#endif
