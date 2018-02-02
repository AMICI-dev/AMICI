#ifndef AMICI_INTERFACE_CPP_H
#define AMICI_INTERFACE_CPP_H

#include "amici_defines.h"
#include "include/amici_misc.h"
#include <include/edata.h>
#include <include/rdata.h>

namespace amici {

/**
 * getSimulationResults runs the forward an backwards simulation and returns
 * results in a ReturnData struct
 *
 * @param[in] edata pointer to the experimental data struct @type ExpData
 * @return rdata data struct with simulation results @type ReturnData
 */

ReturnData *getSimulationResults(Model &model, const ExpData *edata, Solver &solver);

} // namespace amici
#endif
