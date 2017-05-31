#ifndef amici_h
#define amici_h
#include <include/symbolic_functions.h>
#include <include/udata.h>
#include <include/rdata.h>
#include <include/edata.h>
#include <include/tdata.h>
#include <stdbool.h>

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
#define AMI_DATA_RETURN           1
#define AMI_NORMAL                1
#define AMI_ONE_STEP              2

void runAmiciSimulation(UserData *udata, const ExpData *edata, ReturnData *rdata, int *pstatus);

void invalidateReturnData(UserData* udata, ReturnData* rdata);

void *setupAMI(int *status, UserData *udata, TempData *tdata);
void setupAMIB(int *status, void *ami_mem, UserData *udata, TempData *tdata);

void getDataSensisFSA(int *status, int it, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata);
void getDataSensisASA(int *status, int it, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata);

void getEventSensisFSA(int *status, int ie, void *ami_mem, UserData *udata, ReturnData *rdata, TempData *tdata);
void getEventSensisASA(int *status, int ie, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata);

void getEventSigma(int *status, int ie, int iz, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata);
void getEventObjective(int *status, int ie, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata);

void getDataOutput(int *status, int it, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata);
void getEventOutput(int *status, realtype *tlastroot, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata);
void fillEventOutput(int *status, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata);

void handleEvent(int *status, int *iroot, realtype *tlastroot, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata, int seflag);
void handleDataPoint(int *status, int it, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata);
void handleDataPointB(int *status, int it, void *ami_mem, UserData *udata, ReturnData *rdata, TempData *tdata);
void handleEventB(int *status, int iroot, void *ami_mem, UserData *udata, TempData *tdata);

void applyEventBolus(int *status, void *ami_mem, UserData *udata, TempData *tdata);
void applyEventSensiBolusFSA(int *status, void *ami_mem, UserData *udata, TempData *tdata);

realtype getTnext(realtype *troot, int iroot, realtype *tdata, int it, UserData *udata);

void initHeaviside(int *status, UserData *udata, TempData *tdata);
void updateHeaviside(int *status, UserData *udata, TempData *tdata);
void updateHeavisideB(int *status, int iroot, UserData *udata, TempData *tdata);

void getDiagnosis(int *status,int it, void *ami_mem, UserData *udata, ReturnData *rdata);
void getDiagnosisB(int *status,int it, void *ami_mem, UserData *udata, ReturnData *rdata, TempData *tdata);

int workForwardProblem(UserData *udata, TempData *tdata, ReturnData *rdata, const ExpData *edata, int* status, void *ami_mem, int* iroot);
int workBackwardProblem(UserData *udata, TempData *tdata, ReturnData *rdata, const ExpData *edata, int *status, void *ami_mem, int *iroot, booleantype *setupBdone);
void storeJacobianAndDerivativeInReturnData(UserData *udata, TempData *tdata, ReturnData *rdata);
void freeTempDataAmiMem(UserData *udata, TempData *tdata, void *ami_mem, booleantype setupBdone, int status);
void applyChainRuleFactorToSimulationResults(const UserData *udata, ReturnData *rdata, const ExpData *edata);
void unscaleParameters(UserData *udata);

#endif /* amici_h */
