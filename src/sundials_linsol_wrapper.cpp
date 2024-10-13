#include <amici/sundials_linsol_wrapper.h>

#include <amici/exception.h>
#include <sunlinsol/sunlinsol_band.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_klu.h>
#include <sunlinsol/sunlinsol_pcg.h>
#ifdef SUNDIALS_SUPERLUMT
#include <sunlinsol/sunlinsol_superlumt.h>
#endif
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

#include <utility>

namespace amici {

SUNLinSolWrapper::SUNLinSolWrapper(SUNLinearSolver linsol)
    : solver_(linsol) {}

SUNLinSolWrapper::SUNLinSolWrapper(
    SUNLinearSolver linsol, SUNMatrixWrapper const& A
)
    : solver_(linsol)
    , A_(A) {}

SUNLinSolWrapper::~SUNLinSolWrapper() {
    if (solver_)
        SUNLinSolFree(solver_);
}

SUNLinSolWrapper::SUNLinSolWrapper(SUNLinSolWrapper&& other) noexcept {
    std::swap(solver_, other.solver_);
    std::swap(A_, other.A_);
}

SUNLinSolWrapper& SUNLinSolWrapper::operator=(SUNLinSolWrapper&& other
) noexcept {
    std::swap(solver_, other.solver_);
    std::swap(A_, other.A_);
    return *this;
}

SUNLinearSolver SUNLinSolWrapper::get() const { return solver_; }

SUNLinearSolver_Type SUNLinSolWrapper::getType() const {
    return SUNLinSolGetType(solver_);
}

int SUNLinSolWrapper::initialize() {
    auto res = SUNLinSolInitialize(solver_);
    if (res != SUN_SUCCESS)
        throw AmiException(
            "Linear solver initialization failed with code %d", res
        );
    return res;
}

void SUNLinSolWrapper::setup() const {
    auto res = SUNLinSolSetup(solver_, A_.get());
    if (res != SUN_SUCCESS)
        throw AmiException("Linear solver setup failed with code %d", res);
}

int SUNLinSolWrapper::solve(N_Vector x, N_Vector b, realtype tol) const {
    return SUNLinSolSolve(solver_, A_.get(), x, b, tol);
}

long SUNLinSolWrapper::getLastFlag() const {
    return gsl::narrow<long>(SUNLinSolLastFlag(solver_));
}

int SUNLinSolWrapper::space(long* lenrwLS, long* leniwLS) const {
    return SUNLinSolSpace(solver_, lenrwLS, leniwLS);
}

SUNMatrixWrapper& SUNLinSolWrapper::getMatrix() { return A_; }

SUNNonLinSolWrapper::SUNNonLinSolWrapper(SUNNonlinearSolver sol)
    : solver(sol) {}

SUNNonLinSolWrapper::~SUNNonLinSolWrapper() {
    if (solver)
        SUNNonlinSolFree(solver);
}

SUNNonLinSolWrapper::SUNNonLinSolWrapper(SUNNonLinSolWrapper&& other) noexcept {
    std::swap(solver, other.solver);
}

SUNNonLinSolWrapper& SUNNonLinSolWrapper::operator=(SUNNonLinSolWrapper&& other
) noexcept {
    std::swap(solver, other.solver);
    return *this;
}

SUNNonlinearSolver SUNNonLinSolWrapper::get() const { return solver; }

SUNNonlinearSolver_Type SUNNonLinSolWrapper::getType() const {
    return SUNNonlinSolGetType(solver);
}

int SUNNonLinSolWrapper::setup(N_Vector y, void* mem) {
    auto res = SUNNonlinSolSetup(solver, y, mem);
    if (res != SUN_SUCCESS)
        throw AmiException("Nonlinear solver setup failed with code %d", res);
    return res;
}

int SUNNonLinSolWrapper::Solve(
    N_Vector y0, N_Vector y, N_Vector w, realtype tol, bool callLSetup,
    void* mem
) {
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

int SUNNonLinSolWrapper::setConvTestFn(
    SUNNonlinSolConvTestFn CTestFn, void* ctest_data
) {
    return SUNNonlinSolSetConvTestFn(solver, CTestFn, ctest_data);
}

int SUNNonLinSolWrapper::setMaxIters(int maxiters) {
    return SUNNonlinSolSetMaxIters(solver, maxiters);
}

long SUNNonLinSolWrapper::getNumIters() const {
    long int niters = -1;
    auto res = SUNNonlinSolGetNumIters(solver, &niters);
    if (res != SUN_SUCCESS) {
        throw AmiException("SUNNonlinSolGetNumIters failed with code %d", res);
    }
    return niters;
}

int SUNNonLinSolWrapper::getCurIter() const {
    int iter = -1;
    auto res = SUNNonlinSolGetCurIter(solver, &iter);
    if (res != SUN_SUCCESS) {
        throw AmiException("SUNNonlinSolGetCurIter failed with code %d", res);
    }
    return iter;
}

long SUNNonLinSolWrapper::getNumConvFails() const {
    long int nconvfails = -1;
    auto res = SUNNonlinSolGetNumConvFails(solver, &nconvfails);
    if (res != SUN_SUCCESS) {
        throw AmiException(
            "SUNNonlinSolGetNumConvFails failed with code %d", res
        );
    }
    return nconvfails;
}

void SUNNonLinSolWrapper::initialize() {
    int status = SUNNonlinSolInitialize(solver);
    if (status != SUN_SUCCESS)
        throw AmiException(
            "Nonlinear solver initialization failed with code %d", status
        );
}

SUNLinSolBand::SUNLinSolBand(N_Vector x, SUNMatrixWrapper A)
    : SUNLinSolWrapper(SUNLinSol_Band(x, A, x->sunctx)) {
    if (!solver_)
        throw AmiException("Failed to create solver.");
}

SUNLinSolBand::SUNLinSolBand(AmiVector const& x, int ubw, int lbw)
    : SUNLinSolWrapper(
          nullptr, SUNMatrixWrapper(x.getLength(), ubw, lbw, x.get_ctx())
      ) {
    solver_
        = SUNLinSol_Band(const_cast<N_Vector>(x.getNVector()), A_, x.get_ctx());
    if (!solver_)
        throw AmiException("Failed to create solver.");
}

SUNLinSolDense::SUNLinSolDense(AmiVector const& x)
    : SUNLinSolWrapper(
          nullptr, SUNMatrixWrapper(x.getLength(), x.getLength(), x.get_ctx())
      ) {
    solver_ = SUNLinSol_Dense(
        const_cast<N_Vector>(x.getNVector()), A_, x.get_ctx()
    );
    if (!solver_)
        throw AmiException("Failed to create solver.");
}

SUNLinSolKLU::SUNLinSolKLU(N_Vector x, SUNMatrixWrapper A)
    : SUNLinSolWrapper(SUNLinSol_KLU(x, A, x->sunctx)) {
    if (!solver_)
        throw AmiException("Failed to create solver.");
}

SUNLinSolKLU::SUNLinSolKLU(
    AmiVector const& x, int nnz, int sparsetype, StateOrdering ordering
)
    : SUNLinSolWrapper(
          nullptr,
          SUNMatrixWrapper(
              x.getLength(), x.getLength(), nnz, sparsetype, x.get_ctx()
          )
      ) {
    solver_
        = SUNLinSol_KLU(const_cast<N_Vector>(x.getNVector()), A_, x.get_ctx());
    if (!solver_)
        throw AmiException("Failed to create solver.");

    setOrdering(ordering);
}

void SUNLinSolKLU::reInit(int nnz, int reinit_type) {
    int status = SUNLinSol_KLUReInit(solver_, A_, nnz, reinit_type);
    if (status != SUN_SUCCESS)
        throw AmiException("SUNLinSol_KLUReInit failed with %d", status);
}

void SUNLinSolKLU::setOrdering(StateOrdering ordering) {
    auto status = SUNLinSol_KLUSetOrdering(solver_, static_cast<int>(ordering));
    if (status != SUN_SUCCESS)
        throw AmiException("SUNLinSol_KLUSetOrdering failed with %d", status);
}

bool SUNLinSolKLU::is_singular() const {
    // adapted from SUNLinSolSetup_KLU in sunlinsol/klu/sunlinsol_klu.c
    auto content = (SUNLinearSolverContent_KLU)(solver_->content);
    // first cheap check via rcond
    auto status
        = sun_klu_rcond(content->symbolic, content->numeric, &content->common);
    if (status == 0)
        throw AmiException("sun_klu_rcond: %d", content->last_flag);

    auto precision = std::numeric_limits<realtype>::epsilon();

    if (content->common.rcond < precision) {
        // cheap check indicates singular, expensive check via condest
        status = sun_klu_condest(
            SM_INDEXPTRS_S(A_.get()), SM_DATA_S(A_.get()), content->symbolic,
            content->numeric, &content->common
        );
        if (status == 0)
            throw AmiException("sun_klu_condest: %d", content->last_flag);
        return content->common.condest > 1.0 / precision;
    }
    return false;
}

SUNLinSolPCG::SUNLinSolPCG(N_Vector y, int pretype, int maxl)
    : SUNLinSolWrapper(SUNLinSol_PCG(y, pretype, maxl, y->sunctx)) {
    if (!solver_)
        throw AmiException("Failed to create solver.");
}

int SUNLinSolPCG::setATimes(void* A_data, SUNATimesFn ATimes) {
    return SUNLinSolSetATimes_PCG(solver_, A_data, ATimes);
}

int SUNLinSolPCG::setPreconditioner(
    void* P_data, SUNPSetupFn Pset, SUNPSolveFn Psol
) {
    return SUNLinSolSetPreconditioner_PCG(solver_, P_data, Pset, Psol);
}

int SUNLinSolPCG::setScalingVectors(N_Vector s, N_Vector nul) {
    return SUNLinSolSetScalingVectors_PCG(solver_, s, nul);
}

int SUNLinSolPCG::getNumIters() const { return SUNLinSolNumIters_PCG(solver_); }

realtype SUNLinSolPCG::getResNorm() const {
    return SUNLinSolResNorm_PCG(solver_);
}

N_Vector SUNLinSolPCG::getResid() const { return SUNLinSolResid_PCG(solver_); }

SUNLinSolSPBCGS::SUNLinSolSPBCGS(N_Vector x, int pretype, int maxl)
    : SUNLinSolWrapper(SUNLinSol_SPBCGS(x, pretype, maxl, x->sunctx)) {
    if (!solver_)
        throw AmiException("Failed to create solver.");
}

SUNLinSolSPBCGS::SUNLinSolSPBCGS(AmiVector const& x, int pretype, int maxl) {
    solver_ = SUNLinSol_SPBCGS(
        const_cast<N_Vector>(x.getNVector()), pretype, maxl, x.get_ctx()
    );
    if (!solver_)
        throw AmiException("Failed to create solver.");
}

int SUNLinSolSPBCGS::setATimes(void* A_data, SUNATimesFn ATimes) {
    return SUNLinSolSetATimes_SPBCGS(solver_, A_data, ATimes);
}

int SUNLinSolSPBCGS::setPreconditioner(
    void* P_data, SUNPSetupFn Pset, SUNPSolveFn Psol
) {
    return SUNLinSolSetPreconditioner_SPBCGS(solver_, P_data, Pset, Psol);
}

int SUNLinSolSPBCGS::setScalingVectors(N_Vector s, N_Vector nul) {
    return SUNLinSolSetScalingVectors_SPBCGS(solver_, s, nul);
}

int SUNLinSolSPBCGS::getNumIters() const {
    return SUNLinSolNumIters_SPBCGS(solver_);
}

realtype SUNLinSolSPBCGS::getResNorm() const {
    return SUNLinSolResNorm_SPBCGS(solver_);
}

N_Vector SUNLinSolSPBCGS::getResid() const {
    return SUNLinSolResid_SPBCGS(solver_);
}

SUNLinSolSPFGMR::SUNLinSolSPFGMR(AmiVector const& x, int pretype, int maxl)
    : SUNLinSolWrapper(SUNLinSol_SPFGMR(
          const_cast<N_Vector>(x.getNVector()), pretype, maxl, x.get_ctx()
      )) {
    if (!solver_)
        throw AmiException("Failed to create solver.");
}

int SUNLinSolSPFGMR::setATimes(void* A_data, SUNATimesFn ATimes) {
    return SUNLinSolSetATimes_SPFGMR(solver_, A_data, ATimes);
}

int SUNLinSolSPFGMR::setPreconditioner(
    void* P_data, SUNPSetupFn Pset, SUNPSolveFn Psol
) {
    return SUNLinSolSetPreconditioner_SPFGMR(solver_, P_data, Pset, Psol);
}

int SUNLinSolSPFGMR::setScalingVectors(N_Vector s, N_Vector nul) {
    return SUNLinSolSetScalingVectors_SPFGMR(solver_, s, nul);
}

int SUNLinSolSPFGMR::getNumIters() const {
    return SUNLinSolNumIters_SPFGMR(solver_);
}

realtype SUNLinSolSPFGMR::getResNorm() const {
    return SUNLinSolResNorm_SPFGMR(solver_);
}

N_Vector SUNLinSolSPFGMR::getResid() const {
    return SUNLinSolResid_SPFGMR(solver_);
}

SUNLinSolSPGMR::SUNLinSolSPGMR(AmiVector const& x, int pretype, int maxl)
    : SUNLinSolWrapper(SUNLinSol_SPGMR(
          const_cast<N_Vector>(x.getNVector()), pretype, maxl, x.get_ctx()
      )) {
    if (!solver_)
        throw AmiException("Failed to create solver.");
}

int SUNLinSolSPGMR::setATimes(void* A_data, SUNATimesFn ATimes) {
    return SUNLinSolSetATimes_SPGMR(solver_, A_data, ATimes);
}

int SUNLinSolSPGMR::setPreconditioner(
    void* P_data, SUNPSetupFn Pset, SUNPSolveFn Psol
) {
    return SUNLinSolSetPreconditioner_SPGMR(solver_, P_data, Pset, Psol);
}

int SUNLinSolSPGMR::setScalingVectors(N_Vector s, N_Vector nul) {
    return SUNLinSolSetScalingVectors_SPGMR(solver_, s, nul);
}

int SUNLinSolSPGMR::getNumIters() const {
    return SUNLinSolNumIters_SPGMR(solver_);
}

realtype SUNLinSolSPGMR::getResNorm() const {
    return SUNLinSolResNorm_SPGMR(solver_);
}

N_Vector SUNLinSolSPGMR::getResid() const {
    return SUNLinSolResid_SPGMR(solver_);
}

SUNLinSolSPTFQMR::SUNLinSolSPTFQMR(N_Vector x, int pretype, int maxl)
    : SUNLinSolWrapper(SUNLinSol_SPTFQMR(x, pretype, maxl, x->sunctx)) {
    if (!solver_)
        throw AmiException("Failed to create solver.");
}

SUNLinSolSPTFQMR::SUNLinSolSPTFQMR(AmiVector const& x, int pretype, int maxl) {
    solver_ = SUNLinSol_SPTFQMR(
        const_cast<N_Vector>(x.getNVector()), pretype, maxl, x.get_ctx()
    );
    if (!solver_)
        throw AmiException("Failed to create solver.");
}

int SUNLinSolSPTFQMR::setATimes(void* A_data, SUNATimesFn ATimes) {
    return SUNLinSolSetATimes_SPTFQMR(solver_, A_data, ATimes);
}

int SUNLinSolSPTFQMR::setPreconditioner(
    void* P_data, SUNPSetupFn Pset, SUNPSolveFn Psol
) {
    return SUNLinSolSetPreconditioner_SPTFQMR(solver_, P_data, Pset, Psol);
}

int SUNLinSolSPTFQMR::setScalingVectors(N_Vector s, N_Vector nul) {
    return SUNLinSolSetScalingVectors_SPTFQMR(solver_, s, nul);
}

int SUNLinSolSPTFQMR::getNumIters() const {
    return SUNLinSolNumIters_SPTFQMR(solver_);
}

realtype SUNLinSolSPTFQMR::getResNorm() const {
    return SUNLinSolResNorm_SPTFQMR(solver_);
}

N_Vector SUNLinSolSPTFQMR::getResid() const {
    return SUNLinSolResid_SPTFQMR(solver_);
}

SUNNonLinSolNewton::SUNNonLinSolNewton(N_Vector x)
    : SUNNonLinSolWrapper(SUNNonlinSol_Newton(x, x->sunctx)) {}

SUNNonLinSolNewton::SUNNonLinSolNewton(int count, N_Vector x)
    : SUNNonLinSolWrapper(SUNNonlinSol_NewtonSens(count, x, x->sunctx)) {
    if (!solver)
        throw(AmiException("SUNNonlinSol_NewtonSens failed"));
}

int SUNNonLinSolNewton::getSysFn(SUNNonlinSolSysFn* SysFn) const {
    return SUNNonlinSolGetSysFn_Newton(solver, SysFn);
}

SUNNonLinSolFixedPoint::SUNNonLinSolFixedPoint(const_N_Vector x, int m)
    : SUNNonLinSolWrapper(
          SUNNonlinSol_FixedPoint(const_cast<N_Vector>(x), m, x->sunctx)
      ) {}

SUNNonLinSolFixedPoint::SUNNonLinSolFixedPoint(
    int count, const_N_Vector x, int m
)
    : SUNNonLinSolWrapper(SUNNonlinSol_FixedPointSens(
          count, const_cast<N_Vector>(x), m, x->sunctx
      )) {}

int SUNNonLinSolFixedPoint::getSysFn(SUNNonlinSolSysFn* SysFn) const {
    return SUNNonlinSolGetSysFn_FixedPoint(solver, SysFn);
}

#ifdef SUNDIALS_SUPERLUMT

SUNLinSolSuperLUMT::SUNLinSolSuperLUMT(
    N_Vector x, SUNMatrixWrapper A, int numThreads
)
    : SUNLinSolWrapper(SUNLinSol_SuperLUMT(x, A, numThreads), A) {
    if (!solver)
        throw AmiException("Failed to create solver.");
}

SUNLinSolSuperLUMT::SUNLinSolSuperLUMT(
    AmiVector const& x, int nnz, int sparsetype,
    SUNLinSolSuperLUMT::StateOrdering ordering
)
    : SUNLinSolWrapper(
          nullptr,
          SUNMatrixWrapper(x.getLength(), x.getLength(), nnz, sparsetype)
      ) {
    int numThreads = 1;
    if (auto env = std::getenv("AMICI_SUPERLUMT_NUM_THREADS")) {
        numThreads = std::max(1, std::stoi(env));
    }

    solver = SUNLinSol_SuperLUMT(x.getNVector(), A.get(), numThreads);
    if (!solver)
        throw AmiException("Failed to create solver.");

    setOrdering(ordering);
}

SUNLinSolSuperLUMT::SUNLinSolSuperLUMT(
    AmiVector const& x, int nnz, int sparsetype, StateOrdering ordering,
    int numThreads
)
    : SUNLinSolWrapper(
          nullptr,
          SUNMatrixWrapper(x.getLength(), x.getLength(), nnz, sparsetype)
      ) {
    solver = SUNLinSol_SuperLUMT(x.getNVector(), A.get(), numThreads);
    if (!solver)
        throw AmiException("Failed to create solver.");

    setOrdering(ordering);
}

void SUNLinSolSuperLUMT::setOrdering(StateOrdering ordering) {
    auto status
        = SUNLinSol_SuperLUMTSetOrdering(solver, static_cast<int>(ordering));
    if (status != SUNLS_SUCCESS)
        throw AmiException(
            "SUNLinSol_SuperLUMTSetOrdering failed with %d", status
        );
}

#endif

} // namespace amici
