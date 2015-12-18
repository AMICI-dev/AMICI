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

UserData setupUserData(const mxArray *prhs[]);
ReturnData setupReturnData(const mxArray *prhs[], void *user_data);
ExpData setupExpData(const mxArray *prhs[], void *user_data);

void *setupAMI(int *status, void *user_data, void *temp_data);
void setupAMIB(int *status,void *ami_mem, void *user_data, void *temp_data);

void getDataSensisFSA(int *status, int it, void *ami_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data);
void getDataSensisASA(int *status, int it, void *ami_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data);

void getEventSensisFSA(int *status, int ie, void *ami_mem, void  *user_data, void *return_data, void *temp_data);
void getEventSensisASA(int *status, int ie, void *ami_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data);

void getEventSigma(int *status, int ie, int iz, void *ami_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data);
void getEventObjective(int *status, int ie, void *ami_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data);

void getDataOutput(int *status, int it, void *ami_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data);
int getEventOutput(int *status, realtype *tlastroot, void *ami_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data);
void fillEventOutput(int *status, void *ami_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data);

void handleEvent(int *status, int iroot, realtype *tlastroot, void *ami_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data);
void handleDataPoint(int *status, int it, void *ami_mem, void  *user_data, void *return_data, void *exp_data, void *temp_data);
void handleDataPointB(int *status, int it, void *ami_mem, void  *user_data, void *return_data, void *temp_data);
void handleEventB(int *status, int iroot, void *ami_mem, void  *user_data, void *temp_data);

void applyEventBolus(int *status, void *ami_mem, void  *user_data, void *temp_data);
void applyEventSensiBolusFSA(int *status, void *ami_mem, void  *user_data, void *temp_data);

realtype getTnext(realtype troot, int iroot, realtype tdata, int it, void *user_data);

void initHeaviside(int *status, void  *user_data, void *temp_data);
void updateHeaviside(int *status, void  *user_data, void *temp_data);
void updateHeavisideB(int *status, int iroot, void  *user_data, void *temp_data);

void getDiagnosis(int *status,int it, void *ami_mem, void  *user_data, void *return_data);
void getDiagnosisB(int *status,int it, void *ami_mem, void  *user_data, void *return_data, void *temp_data);



#endif /* amici_symbolic_functions_h */