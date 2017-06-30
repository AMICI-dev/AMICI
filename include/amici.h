#ifndef amici_h
#define amici_h
#include <include/symbolic_functions.h>
#include <include/udata.h>
#include <include/rdata.h>
#include <include/edata.h>
#include <include/tdata.h>
#include <cstdbool>

/* linear solvers */
#define AMICI_DENSE       1
#define AMICI_BAND        2
#define AMICI_LAPACKDENSE 3
#define AMICI_LAPACKBAND  4
#define AMICI_DIAG        5
#define AMICI_SPGMR       6
#define AMICI_SPBCG       7
#define AMICI_SPTFQMR     8
#define AMICI_KLU         9

#define AMICI_NORMAL      1
#define AMICI_LOGNORMAL   2
#define AMICI_ONEOUTPUT   5


#define AMICI_ERROR_UDATA         -99
#define AMICI_ERROR_EDATA         -98
#define AMICI_ERROR_RDATA         -97
#define AMICI_ERROR_TDATA         -96
#define AMICI_ERROR_SETUP         -95
#define AMICI_ERROR_SETUPB        -94
#define AMICI_ERROR_NOTHINGTODO   -93
#define AMICI_ERROR_FSA           -92
#define AMICI_ERROR_ASA           -91
#define AMICI_ERROR_SA            -90
#define AMICI_ERROR_SS            -89
#define AMICI_ERROR_DATA          -88
#define AMICI_ERROR_EVENT         -87
#define AMICI_SUCCESS               0
#define AMICI_ROOT_RETURN           2
#define AMICI_DATA_RETURN           1

#define AMICI_NORMAL                1
#define AMICI_ONE_STEP              2

typedef enum {AMICI_BLAS_RowMajor=101, AMICI_BLAS_ColMajor=102} AMICI_BLAS_LAYOUT;
typedef enum {AMICI_BLAS_NoTrans=111, AMICI_BLAS_Trans=112, AMICI_BLAS_ConjTrans=113} AMICI_BLAS_TRANSPOSE;

void errMsgIdAndTxt(
                            const char * identifier, /* string with error message identifier */
                            const char * err_msg,    /* string with error message printf-style format */
                            ...                      /* any additional arguments */
);

void warnMsgIdAndTxt(
                             const char * identifier, /* string with error message identifier */
                             const char * err_msg,    /* string with error message printf-style format */
                             ...                      /* any additional arguments */
);

int runAmiciSimulation(UserData *udata, const ExpData *edata, ReturnData *rdata);

void invalidateReturnData(UserData* udata, ReturnData* rdata);

void *setupAMI(UserData *udata, TempData *tdata);
int setupAMIB(void *ami_mem, UserData *udata, TempData *tdata);

int getDataSensisFSA(int it, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata);
int getDataSensisASA(int it, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata);

int getEventSensisFSA(int ie, void *ami_mem, UserData *udata, ReturnData *rdata, TempData *tdata);
int getEventSensisASA(int ie, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata);

int getEventSigma(int ie, int iz, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata);
int getEventObjective(int ie, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata);

int getDataOutput(int it, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata);
int getEventOutput(realtype *tlastroot, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata);
int fillEventOutput(void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata);

int handleEvent(int *iroot, realtype *tlastroot, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata, int seflag);
int handleDataPoint(int it, void *ami_mem, UserData *udata, ReturnData *rdata, const ExpData *edata, TempData *tdata);
int handleDataPointB(int it, void *ami_mem, UserData *udata, ReturnData *rdata, TempData *tdata);
int handleEventB(int iroot, void *ami_mem, UserData *udata, TempData *tdata);

int applyEventBolus(void *ami_mem, UserData *udata, TempData *tdata);
int applyEventSensiBolusFSA(void *ami_mem, UserData *udata, TempData *tdata);

realtype getTnext(realtype *troot, int iroot, realtype *tdata, int it, UserData *udata);

int initHeaviside(UserData *udata, TempData *tdata);
int updateHeaviside(UserData *udata, TempData *tdata);
int updateHeavisideB(int iroot, UserData *udata, TempData *tdata);

int getDiagnosis(int it, void *ami_mem, UserData *udata, ReturnData *rdata);
int getDiagnosisB(int it, void *ami_mem, UserData *udata, ReturnData *rdata, TempData *tdata);

int workForwardProblem(UserData *udata, TempData *tdata, ReturnData *rdata, const ExpData *edata, void *ami_mem, int* iroot);
int workBackwardProblem(UserData *udata, TempData *tdata, ReturnData *rdata, const ExpData *edata, void *ami_mem, int *iroot);
int storeJacobianAndDerivativeInReturnData(UserData *udata, TempData *tdata, ReturnData *rdata);
int applyChainRuleFactorToSimulationResults(const UserData *udata, ReturnData *rdata, const ExpData *edata);
int unscaleParameters(UserData *udata);

void amici_dgemv(AMICI_BLAS_LAYOUT layout,
                 AMICI_BLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY);

void amici_dgemm(AMICI_BLAS_LAYOUT layout, AMICI_BLAS_TRANSPOSE TransA,
                 AMICI_BLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);

#endif /* amici_h */
