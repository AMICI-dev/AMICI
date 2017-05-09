#ifndef AMICI_INTERFACE_MATLAB_H
#define AMICI_INTERFACE_MATLAB_H

#include <mex.h>
#include "include/amici.h"

/**
 * @brief userDataFromMatlabCall extracts information from the matlab call and returns the corresponding UserData struct
 * @param[in] prhs: pointer to the array of input arguments @type mxArray
 * @return udata: struct containing all provided user data @type UserData
 */
UserData *userDataFromMatlabCall(const mxArray *prhs[]);

/**
 * setupReturnData initialises the return data struct
 * @param[in] plhs user input @type mxArray
 * @param[in] udata pointer to the user data struct @type UserData
 * @param[out] pstatus pointer to the flag indicating the execution status @type double
 * @return rdata: return data struct @type ReturnData
 */
ReturnData *setupReturnData(mxArray *plhs[], const UserData *udata, double *pstatus);

/**
 * setupExpData initialises the experimental data struct
 * @param[in] prhs user input @type *mxArray
 * @param[in] udata pointer to the user data struct @type UserData
 * @return edata: experimental data struct @type ExpData
 */
ExpData *setupExpData(const mxArray *prhs[], UserData *udata, int *status);


#endif
