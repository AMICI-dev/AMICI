#ifndef AMICI_INTERFACE_CPP_H
#define AMICI_INTERFACE_CPP_H

#include "include/amici_misc.h"

#include <include/udata.h>
#include <include/edata.h>
#include <include/rdata.h>

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

/**
 * getSimulationResults runs the forward an backwards simulation and returns results in a ReturnData struct
 *
 * @param[in] udata pointer to the user data struct @type UserData
 * @param[in] edata pointer to the experimental data struct @type ExpData
 * @param[out] pstatus flag indicating success of execution @type *int
 * @return rdata data struct with simulation results @type ReturnData
 */

EXTERNC ReturnData *getSimulationResults(UserData *udata, const ExpData *edata);

#endif
