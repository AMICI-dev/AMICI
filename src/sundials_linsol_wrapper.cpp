#include <amici/sundials_linsol_wrapper.h>
#include <amici/exception.h>

#include <utility>
#include <new> // bad_alloc

namespace amici {

SUNLinSolWrapper::SUNLinSolWrapper(SUNLinearSolver linsol_)
    : linsol(linsol_)
{
}

SUNLinSolWrapper::~SUNLinSolWrapper()
{
    if(linsol)
        SUNLinSolFree(linsol);
}

SUNLinSolWrapper::SUNLinSolWrapper(SUNLinSolWrapper &&other) noexcept
{
    std::swap(linsol, other.linsol);
}

SUNLinearSolver SUNLinSolWrapper::get() const
{
    return linsol;
}

SUNLinearSolver_Type SUNLinSolWrapper::getType() const
{
    return SUNLinSolGetType(linsol);
}

int SUNLinSolWrapper::initialize()
{
    auto res = SUNLinSolInitialize(linsol);
    if(res != SUNLS_SUCCESS)
        throw AmiException("Solver initialization failed with code %d", res);
    return res;
}

int SUNLinSolWrapper::setup(SUNMatrix A)
{
    auto res = SUNLinSolSetup(linsol, A);
    if(res != SUNLS_SUCCESS)
        throw AmiException("Solver setup failed with code %d", res);
    return res;
}

int SUNLinSolWrapper::setup(SUNMatrixWrapper A)
{
    return setup(A.get());
}

int SUNLinSolWrapper::Solve(SUNMatrix A, N_Vector x, N_Vector b, realtype tol)
{
    // TODO: tol as member?
    return SUNLinSolSolve(linsol, A, x, b, tol);
}

long SUNLinSolWrapper::getLastFlag()
{
    return SUNLinSolLastFlag(linsol);
}

int SUNLinSolWrapper::space(long *lenrwLS, long *leniwLS)
{
    return SUNLinSolSpace(linsol, lenrwLS, leniwLS);
}

SUNMatrix SUNLinSolWrapper::getMatrix() const
{
    return nullptr;
}

SUNNonLinSolWrapper::SUNNonLinSolWrapper(SUNNonlinearSolver sol_)
    :solver(sol_)
{
}

SUNNonLinSolWrapper::~SUNNonLinSolWrapper()
{
    if(solver)
        SUNNonlinSolFree(solver);
}


SUNNonLinSolWrapper::SUNNonLinSolWrapper(SUNNonLinSolWrapper &&other) noexcept
{
    std::swap(solver, other.solver);
}


SUNNonLinSolWrapper &SUNNonLinSolWrapper::operator=(SUNNonLinSolWrapper &&other) noexcept
{
    std::swap(solver, other.solver);
    return *this;

}

SUNNonlinearSolver SUNNonLinSolWrapper::get() const
{
    return solver;
}

SUNNonlinearSolver_Type SUNNonLinSolWrapper::getType() const
{
    return SUNNonlinSolGetType(solver);
}

int SUNNonLinSolWrapper::setup(N_Vector y, void *mem)
{
    auto res = SUNNonlinSolSetup(solver, y, mem);
    if(res != SUN_NLS_SUCCESS)
        throw AmiException("Nonlinear solver setup failed with code %d", res);
    return res;
}

int SUNNonLinSolWrapper::Solve(N_Vector y0, N_Vector y, N_Vector w, realtype tol, int callLSetup, void *mem)
{
    return SUNNonlinSolSolve(solver, y0, y, w, tol, callLSetup, mem);
}

int SUNNonLinSolWrapper::setSysFn(SUNNonlinSolSysFn SysFn) {
    return SUNNonlinSolSetSysFn(solver, SysFn);
}

int SUNNonLinSolWrapper::setLSetupFn(SUNNonlinSolLSetupFn SetupFn) {
    return SUNNonlinSolSetLSetupFn(solver, SetupFn);
}

int SUNNonLinSolWrapper::setLSolveFn(SUNNonlinSolLSolveFn SolveFn) {
    return SUNNonlinSolSetLSolveFn(solver, SolveFn);
}

int SUNNonLinSolWrapper::setConvTestFn(SUNNonlinSolConvTestFn CTestFn) {
    return SUNNonlinSolSetConvTestFn(solver, CTestFn);
}

int SUNNonLinSolWrapper::setMaxIters(int maxiters) {
    return SUNNonlinSolSetMaxIters(solver, maxiters);
}

long SUNNonLinSolWrapper::getNumIters() {
    long int niters = -1;
    auto res = SUNNonlinSolGetNumIters(solver, &niters);
    if(res != SUN_NLS_SUCCESS) {
        throw AmiException("SUNNonlinSolGetNumIters failed with code %d", res);
    }
    return niters;
}

int SUNNonLinSolWrapper::getCurIter() {
    int iter = -1;
    auto res = SUNNonlinSolGetCurIter(solver, &iter);
    if(res != SUN_NLS_SUCCESS) {
        throw AmiException("SUNNonlinSolGetCurIter failed with code %d", res);
    }
    return iter;
}

long SUNNonLinSolWrapper::getNumConvFails() {
    long int nconvfails = -1;
    auto res = SUNNonlinSolGetNumConvFails(solver, &nconvfails);
    if(res != SUN_NLS_SUCCESS) {
        throw AmiException("SUNNonlinSolGetNumConvFails failed with code %d", res);
    }
    return nconvfails;
}

void SUNNonLinSolWrapper::initialize()
{
    int status = SUNNonlinSolInitialize(solver);
    if(status != SUN_NLS_SUCCESS)
        throw AmiException("Nonlinear solver initialization failed with code %d", status);
}

SUNLinSolBand::SUNLinSolBand(N_Vector y, SUNMatrix A)
    : SUNLinSolWrapper(SUNLinSol_Band(y, A))
{
    if(!linsol)
        throw AmiException("Failed to create solver.");
    initialize();
}

SUNLinSolBand::SUNLinSolBand(const AmiVector &x, int ubw, int lbw)
{
    A = SUNMatrixWrapper(x.getLength(), ubw, lbw);
    linsol = SUNLinSol_Band(x.getNVector(), A.get());
}

SUNMatrix SUNLinSolBand::getMatrix() const
{
    return A.get();
}

SUNLinSolDense::SUNLinSolDense(const AmiVector &x)
{
    A = SUNMatrixWrapper(x.getLength(), x.getLength());
    linsol = SUNLinSol_Dense(x.getNVector(), A.get());
    if(!linsol)
        throw AmiException("Failed to create solver.");
    initialize();
}

SUNMatrix SUNLinSolDense::getMatrix() const
{
    return A.get();
}

SUNLinSolKLU::SUNLinSolKLU(N_Vector y, SUNMatrix A)
    : SUNLinSolWrapper(SUNLinSol_KLU(y, A))
{
    if(!linsol)
        throw AmiException("Failed to create solver.");
    initialize();
}

SUNLinSolKLU::SUNLinSolKLU(const AmiVector &x, int nnz, int sparsetype, StateOrdering ordering)
{
    A = SUNMatrixWrapper(x.getLength(), x.getLength(), nnz, sparsetype);
    linsol = SUNLinSol_KLU(x.getNVector(), A.get());

    if(!linsol)
        throw AmiException("Failed to create solver.");
    initialize();

    setOrdering(ordering);
}

SUNMatrix SUNLinSolKLU::getMatrix() const
{
    return A.get();
}

void SUNLinSolKLU::reInit(int nnz, int reinit_type) {
    int status = SUNLinSol_KLUReInit(linsol, A.get(), nnz, reinit_type);
    if(status != SUNLS_SUCCESS)
        throw AmiException("SUNLinSol_KLUReInit failed with %d", status);
}

void SUNLinSolKLU::setOrdering(StateOrdering ordering) {
    auto status = SUNLinSol_KLUSetOrdering(linsol, static_cast<int>(ordering));
    if(status != SUNLS_SUCCESS)
        throw AmiException("SUNLinSol_KLUSetOrdering failed with %d", status);

}

SUNLinSolPCG::SUNLinSolPCG(N_Vector y, int pretype, int maxl)
    : SUNLinSolWrapper(SUNLinSol_PCG(y, pretype, maxl))
{
    if(!linsol)
        throw AmiException("Failed to create solver.");
    initialize();
}

int SUNLinSolPCG::setATimes(void *A_data, ATimesFn ATimes)
{
    return SUNLinSolSetATimes_PCG(linsol, A_data, ATimes);
}

int SUNLinSolPCG::setPreconditioner(void *P_data, PSetupFn Pset, PSolveFn Psol) {
    return SUNLinSolSetPreconditioner_PCG(linsol, P_data, Pset, Psol);
}

int SUNLinSolPCG::setScalingVectors(N_Vector s, N_Vector nul) {
    return SUNLinSolSetScalingVectors_PCG(linsol, s, nul);
}

int SUNLinSolPCG::getNumIters() {
    return SUNLinSolNumIters_PCG(linsol);
}

realtype SUNLinSolPCG::getResNorm() {
    return SUNLinSolResNorm_PCG(linsol);
}

N_Vector SUNLinSolPCG::getResid() {
    return SUNLinSolResid_PCG(linsol);
}

SUNLinSolSPBCGS::SUNLinSolSPBCGS(N_Vector y, int pretype, int maxl)
    : SUNLinSolWrapper(SUNLinSol_SPBCGS(y, pretype, maxl))
{
    if(!linsol)
        throw AmiException("Failed to create solver.");
    initialize();
}

SUNLinSolSPBCGS::SUNLinSolSPBCGS(const AmiVector &x, int pretype, int maxl)
{
    linsol = SUNLinSol_SPBCGS(x.getNVector(), pretype, maxl);
    if(!linsol)
        throw AmiException("Failed to create solver.");
    initialize();
}

int SUNLinSolSPBCGS::setATimes(void *A_data, ATimesFn ATimes)
{
    return SUNLinSolSetATimes_SPBCGS(linsol, A_data, ATimes);
}

int SUNLinSolSPBCGS::setPreconditioner(void *P_data, PSetupFn Pset, PSolveFn Psol) {
    return SUNLinSolSetPreconditioner_SPBCGS(linsol, P_data, Pset, Psol);
}

int SUNLinSolSPBCGS::setScalingVectors(N_Vector s, N_Vector nul) {
    return SUNLinSolSetScalingVectors_SPBCGS(linsol, s, nul);
}

int SUNLinSolSPBCGS::getNumIters() {
    return SUNLinSolNumIters_SPBCGS(linsol);
}

realtype SUNLinSolSPBCGS::getResNorm() {
    return SUNLinSolResNorm_SPBCGS(linsol);
}

N_Vector SUNLinSolSPBCGS::getResid() {
    return SUNLinSolResid_SPBCGS(linsol);
}

SUNLinSolSPFGMR::SUNLinSolSPFGMR(const AmiVector &x, int pretype, int maxl)
    : SUNLinSolWrapper(SUNLinSol_SPFGMR(x.getNVector(), pretype, maxl))
{
    if(!linsol)
        throw AmiException("Failed to create solver.");
    initialize();
}

int SUNLinSolSPFGMR::setATimes(void *A_data, ATimesFn ATimes)
{
    return SUNLinSolSetATimes_SPFGMR(linsol, A_data, ATimes);
}

int SUNLinSolSPFGMR::setPreconditioner(void *P_data, PSetupFn Pset, PSolveFn Psol) {
    return SUNLinSolSetPreconditioner_SPFGMR(linsol, P_data, Pset, Psol);
}

int SUNLinSolSPFGMR::setScalingVectors(N_Vector s, N_Vector nul) {
    return SUNLinSolSetScalingVectors_SPFGMR(linsol, s, nul);
}

int SUNLinSolSPFGMR::getNumIters() {
    return SUNLinSolNumIters_SPFGMR(linsol);
}

realtype SUNLinSolSPFGMR::getResNorm() {
    return SUNLinSolResNorm_SPFGMR(linsol);
}

N_Vector SUNLinSolSPFGMR::getResid() {
    return SUNLinSolResid_SPFGMR(linsol);
}

SUNLinSolSPGMR::SUNLinSolSPGMR(const AmiVector &x, int pretype, int maxl)
    : SUNLinSolWrapper(SUNLinSol_SPGMR(x.getNVector(), pretype, maxl))
{
    if(!linsol)
        throw AmiException("Failed to create solver.");
    initialize();
}

int SUNLinSolSPGMR::setATimes(void *A_data, ATimesFn ATimes)
{
    return SUNLinSolSetATimes_SPGMR(linsol, A_data, ATimes);
}

int SUNLinSolSPGMR::setPreconditioner(void *P_data, PSetupFn Pset, PSolveFn Psol) {
    return SUNLinSolSetPreconditioner_SPGMR(linsol, P_data, Pset, Psol);
}

int SUNLinSolSPGMR::setScalingVectors(N_Vector s, N_Vector nul) {
    return SUNLinSolSetScalingVectors_SPGMR(linsol, s, nul);
}

int SUNLinSolSPGMR::getNumIters() {
    return SUNLinSolNumIters_SPGMR(linsol);
}

realtype SUNLinSolSPGMR::getResNorm() {
    return SUNLinSolResNorm_SPGMR(linsol);
}

N_Vector SUNLinSolSPGMR::getResid() {
    return SUNLinSolResid_SPGMR(linsol);
}

SUNLinSolSPTFQMR::SUNLinSolSPTFQMR(N_Vector y, int pretype, int maxl)
    : SUNLinSolWrapper(SUNLinSol_SPTFQMR(y, pretype, maxl))
{
    if(!linsol)
        throw AmiException("Failed to create solver.");
    initialize();
}

SUNLinSolSPTFQMR::SUNLinSolSPTFQMR(const AmiVector &x, int pretype, int maxl)
{
    linsol = SUNLinSol_SPTFQMR(x.getNVector(), pretype, maxl);
    if(!linsol)
        throw AmiException("Failed to create solver.");
    initialize();
}

int SUNLinSolSPTFQMR::setATimes(void *A_data, ATimesFn ATimes)
{
    return SUNLinSolSetATimes_SPTFQMR(linsol, A_data, ATimes);
}

int SUNLinSolSPTFQMR::setPreconditioner(void *P_data, PSetupFn Pset, PSolveFn Psol) {
    return SUNLinSolSetPreconditioner_SPTFQMR(linsol, P_data, Pset, Psol);
}

int SUNLinSolSPTFQMR::setScalingVectors(N_Vector s, N_Vector nul) {
    return SUNLinSolSetScalingVectors_SPTFQMR(linsol, s, nul);
}

int SUNLinSolSPTFQMR::getNumIters() {
    return SUNLinSolNumIters_SPTFQMR(linsol);
}

realtype SUNLinSolSPTFQMR::getResNorm() {
    return SUNLinSolResNorm_SPTFQMR(linsol);
}

N_Vector SUNLinSolSPTFQMR::getResid() {
    return SUNLinSolResid_SPTFQMR(linsol);
}

SUNNonLinSolNewton::SUNNonLinSolNewton(N_Vector x)
    :SUNNonLinSolWrapper(SUNNonlinSol_Newton(x))
{
    //initialize();
}

SUNNonLinSolNewton::SUNNonLinSolNewton(int count, N_Vector x)
    :SUNNonLinSolWrapper(SUNNonlinSol_NewtonSens(count, x))
{
    if(!solver)
        throw(AmiException("SUNNonlinSol_NewtonSens failed"));
    //initialize();
}

int SUNNonLinSolNewton::getSysFn(SUNNonlinSolSysFn *SysFn) {
    return SUNNonlinSolGetSysFn_Newton(solver, SysFn);
}

SUNNonLinSolFixedPoint::SUNNonLinSolFixedPoint(N_Vector x, int m)
    :SUNNonLinSolWrapper(SUNNonlinSol_FixedPoint(x, m))
{
    //initialize();
}

SUNNonLinSolFixedPoint::SUNNonLinSolFixedPoint(int count, N_Vector x, int m)
    :SUNNonLinSolWrapper(SUNNonlinSol_FixedPointSens(count, x, m))
{
    //initialize();
}

int SUNNonLinSolFixedPoint::getSysFn(SUNNonlinSolSysFn *SysFn) {
    return SUNNonlinSolGetSysFn_FixedPoint(solver, SysFn);
}



} // namespace amici
