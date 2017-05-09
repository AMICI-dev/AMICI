#ifndef AMICI_INTERFACE_MATLAB_H
#define AMICI_INTERFACE_MATLAB_H

#include <mex.h>
#include "include/amici.h"

UserData *setupUserData(const mxArray *prhs[]);
ReturnData *setupReturnData(mxArray *plhs[], UserData *udata, double *pstatus);
ExpData *setupExpData(const mxArray *prhs[], UserData *udata, int *status);


#endif
