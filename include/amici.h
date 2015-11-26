#ifndef amici_h
#define amici_h
#include <include/udata.h>
#include <include/rdata.h>
#include <include/edata.h>
#include <include/tdata.h>

/* sensitivity method */
#define AMI_FSA 1
#define AMI_ASA 2
#define AMI_SS  3

/* linear solvers */
#define AMI_DENSE       1
#define AMI_BAND        2
#define AMI_LAPACKDENSE 3
#define AMI_LAPACKBAND  4
#define AMI_DIAG        5
#define AMI_SPGMR       6
#define AMI_SPBCG       7
#define AMI_SPTFQMR     8
#define AMI_KLU         9

#define AMI_NORMAL      1
#define AMI_LOGNORMAL   2
#define AMI_ONEOUTPUT   5

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