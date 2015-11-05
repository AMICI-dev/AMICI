#ifndef amici_h
#define amici_h
#include <include/udata.h>
#include <include/rdata.h>
#include <include/edata.h>
#include <include/tdata.h>

/* sensitivity method */
#define CW_FSA 1
#define CW_ASA 2

/* linear solvers */
#define CW_DENSE       1
#define CW_BAND        2
#define CW_LAPACKDENSE 3
#define CW_LAPACKBAND  4
#define CW_DIAG        5
#define CW_SPGMR       6
#define CW_SPBCG       7
#define CW_SPTFQMR     8
#define CW_KLU         9

#define LW_NORMAL      1
#define LW_LOGNORMAL   2
#define LW_ONEOUTPUT   5

#define AMI_SUCCESS               0
#define AMI_ROOT_RETURN           2
#define AMI_NORMAL                1
#define AMI_ONE_STEP              2

#include <src/amici.c>

UserData setupUserData(const mxArray *prhs[]);
ReturnData setupReturnData(const mxArray *prhs[], void *user_data);
ExpData setupExpData(const mxArray *prhs[], void *user_data);

void *setupCVode(int *status, void *user_data, void *temp_data);
void setupCVodeB(int *status, void *cvode_mem, void *user_data, void *temp_data);

void getRootDataFSA(int *status, int *nroots, void *cvode_mem, void  *user_data, void *return_data, void *temp_data);
void getRootDataASA(int *status, int *nroots, int *idisc, void *cvode_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data);

void getDiagnosis(int *status, int it, void *cvode_mem, void  *user_data, void *return_data);


#endif /* amici_symbolic_functions_h */