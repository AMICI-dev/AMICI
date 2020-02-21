#include "amici/solver.h"

#include "amici/exception.h"
#include "amici/misc.h"
#include "amici/model.h"
#include "amici/rdata.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <memory>

namespace amici {

Solver::Solver(AmiciApplication *app) : app(app)
{

}

Solver::Solver(const Solver &other)
    : ism(other.ism), lmm(other.lmm), iter(other.iter),
      interpType(other.interpType), maxsteps(other.maxsteps), t(nan("")),
      sensi_meth(other.sensi_meth), stldet(other.stldet),
      ordering(other.ordering), newton_maxsteps(other.newton_maxsteps),
      newton_maxlinsteps(other.newton_maxlinsteps),
      newton_damping_factor_mode(other.newton_damping_factor_mode),
      newton_damping_factor_lower_bound(other.newton_damping_factor_lower_bound),
      requires_preequilibration(other.requires_preequilibration), linsol(other.linsol),
      atol(other.atol), rtol(other.rtol), atol_fsa(other.atol_fsa),
      rtol_fsa(other.rtol_fsa), atolB(other.atolB), rtolB(other.rtolB),
      quad_atol(other.quad_atol), quad_rtol(other.quad_rtol),
      ss_atol(other.ss_atol), ss_rtol(other.ss_rtol),
      ss_atol_sensi(other.ss_atol_sensi), ss_rtol_sensi(other.ss_rtol_sensi),
      maxstepsB(other.maxstepsB), sensi(other.sensi)
{}

int Solver::run(const realtype tout) const {
    setStopTime(tout);
    clock_t starttime = clock();
    int status;
    if (getAdjInitDone()) {
        status = solveF(tout, AMICI_NORMAL, &ncheckPtr);
    } else {
        status = solve(tout, AMICI_NORMAL);
    }
    cpu_time += (realtype)((clock() - starttime) * 1000) / CLOCKS_PER_SEC;
    return status;
}

int Solver::step(const realtype tout) const {
    int status;
    if (getAdjInitDone()) {
        status = solveF(tout, AMICI_ONE_STEP, &ncheckPtr);
    } else {
        status = solve(tout, AMICI_ONE_STEP);
    }
    return status;
}

void Solver::runB(const realtype tout) const {
    clock_t starttime = clock();
    solveB(tout, AMICI_NORMAL);
    cpu_timeB += (realtype)((clock() - starttime) * 1000) / CLOCKS_PER_SEC;
    t = tout;
}

void Solver::setup(const realtype t0, Model *model, const AmiVector &x0,
                   const AmiVector &dx0, const AmiVectorArray &sx0,
                   const AmiVectorArray &sdx0) const {
    if (nx() != model->nx_solver || nplist() != model->nplist() ||
        nquad() != model->nJ * model->nplist()) {
        resetMutableMemory(model->nx_solver, model->nplist(),
                           model->nJ * model->nplist());
    }
    /* Create solver memory object if necessary */
    allocateSolver();
    if (!solverMemory)
        throw AmiException("Failed to allocated solver memory!");

    /* Initialize CVodes/IDAs solver*/
    init(t0, x0, dx0);

    applyTolerances();

    /* Set optional inputs */
    setErrHandlerFn();
    /* attaches userdata*/
    setUserData(model);
    /* specify maximal number of steps */
    setMaxNumSteps(maxsteps);
    /* activates stability limit detection */
    setStabLimDet(stldet);

    rootInit(model->ne);

    initializeLinearSolver(model);
    initializeNonLinearSolver();

    if (sensi >= SensitivityOrder::first && model->nx_solver > 0) {
        auto plist = model->getParameterList();
        if (sensi_meth == SensitivityMethod::forward && !plist.empty()) {
            /* Set sensitivity analysis optional inputs */
            auto par = model->getUnscaledParameters();

            /* Activate sensitivity calculations */
            sensInit1(sx0, sdx0);
            initializeNonLinearSolverSens(model);
            setSensParams(par.data(), nullptr, plist.data());

            applyTolerancesFSA();
        } else if (sensi_meth == SensitivityMethod::adjoint) {
            /* Allocate space for the adjoint computation */
            adjInit();
        }
    }

    setId(model);
    setSuppressAlg(true);
    /* calculate consistent DAE initial conditions (no effect for ODE) */
    if (model->nt() > 1)
        calcIC(model->getTimepoint(1));
}

void Solver::setupB(int *which, const realtype tf, Model *model,
                    const AmiVector &xB0, const AmiVector &dxB0,
                    const AmiVector &xQB0) const {
    if (!solverMemory)
        throw AmiException(
            "Solver for the forward problem must be setup first");

    /* allocate memory for the backward problem */
    allocateSolverB(which);

    /* initialise states */
    binit(*which, tf, xB0, dxB0);

    /* Attach user data */
    setUserDataB(*which, model);

    /* Number of maximal internal steps */
    setMaxNumStepsB(*which, (maxstepsB == 0) ? maxsteps * 100 : maxstepsB);

    initializeLinearSolverB(model, *which);
    initializeNonLinearSolverB(*which);

    /* Initialise quadrature calculation */
    qbinit(*which, xQB0);

    applyTolerancesASA(*which);
    applyQuadTolerancesASA(*which);

    setStabLimDetB(*which, stldet);
}

void Solver::getDiagnosis(const int it, ReturnData *rdata) const {
    long int number;

    if (solverWasCalledF && solverMemory) {
        getNumSteps(solverMemory.get(), &number);
        rdata->numsteps[it] = number;

        getNumRhsEvals(solverMemory.get(), &number);
        rdata->numrhsevals[it] = number;

        getNumErrTestFails(solverMemory.get(), &number);
        rdata->numerrtestfails[it] = number;

        getNumNonlinSolvConvFails(solverMemory.get(), &number);
        rdata->numnonlinsolvconvfails[it] = number;

        getLastOrder(solverMemory.get(), &rdata->order[it]);
    }
}

void Solver::getDiagnosisB(const int it, ReturnData *rdata,
                           const int which) const {
    long int number;

    if (solverWasCalledB && solverMemoryB.at(which)) {
        getNumSteps(solverMemoryB.at(which).get(), &number);
        rdata->numstepsB[it] = number;

        getNumRhsEvals(solverMemoryB.at(which).get(), &number);
        rdata->numrhsevalsB[it] = number;

        getNumErrTestFails(solverMemoryB.at(which).get(), &number);
        rdata->numerrtestfailsB[it] = number;

        getNumNonlinSolvConvFails(solverMemoryB.at(which).get(), &number);
        rdata->numnonlinsolvconvfailsB[it] = number;
    }
}

void Solver::initializeLinearSolver(const Model *model) const {
    switch (linsol) {

        /* DIRECT SOLVERS */

    case LinearSolver::dense:
        linearSolver = std::make_unique<SUNLinSolDense>(x);
        setLinearSolver();
        setDenseJacFn();
        break;

    case LinearSolver::band:
        linearSolver =
            std::make_unique<SUNLinSolBand>(x, model->ubw, model->lbw);
        setLinearSolver();
        setBandJacFn();
        break;

    case LinearSolver::LAPACKDense:
        throw AmiException("Solver currently not supported!");

    case LinearSolver::LAPACKBand:
        throw AmiException("Solver currently not supported!");

    case LinearSolver::diag:
        diag();
        setDenseJacFn();
        break;

        /* ITERATIVE SOLVERS */

    case LinearSolver::SPGMR:
        linearSolver = std::make_unique<SUNLinSolSPGMR>(x);
        setLinearSolver();
        setJacTimesVecFn();
        break;

    case LinearSolver::SPBCG:
        linearSolver = std::make_unique<SUNLinSolSPBCGS>(x);
        setLinearSolver();
        setJacTimesVecFn();
        break;

    case LinearSolver::SPTFQMR:
        linearSolver = std::make_unique<SUNLinSolSPTFQMR>(x);
        setLinearSolver();
        setJacTimesVecFn();
        break;

        /* SPARSE SOLVERS */

    case LinearSolver::KLU:
        linearSolver = std::make_unique<SUNLinSolKLU>(
            x, model->nnz, CSC_MAT,
            static_cast<SUNLinSolKLU::StateOrdering>(getStateOrdering()));
        setLinearSolver();
        setSparseJacFn();
        break;

#ifdef SUNDIALS_SUPERLUMT
    case LinearSolver::SuperLUMT:
        // TODO state ordering
        linearSolver = std::make_unique<SUNLinSolSuperLUMT>(
            *x, model->nnz, CSC_MAT,
            static_cast<SUNLinSolSuperLUMT::StateOrdering>(getStateOrdering()));

        setLinearSolver();
        setSparseJacFn();
        break;
#endif
    default:
        throw AmiException("Invalid choice of solver: %d",
                           static_cast<int>(linsol));
    }
}

void Solver::initializeNonLinearSolver() const {
    switch (iter) {
    case NonlinearSolverIteration::newton:
        nonLinearSolver = std::make_unique<SUNNonLinSolNewton>(x.getNVector());
        break;
    case NonlinearSolverIteration::fixedpoint:
        nonLinearSolver =
            std::make_unique<SUNNonLinSolFixedPoint>(x.getNVector());
        break;
    default:
        throw AmiException("Invalid non-linear solver specified (%d).",
                           static_cast<int>(iter));
    }

    setNonLinearSolver();
}

void Solver::initializeLinearSolverB(const Model *model,
                                     const int which) const {
    switch (linsol) {
    /* DIRECT SOLVERS */
    case LinearSolver::dense:
        linearSolverB = std::make_unique<SUNLinSolDense>(xB);
        setLinearSolverB(which);
        setDenseJacFnB(which);
        break;

    case LinearSolver::band:
        linearSolverB =
            std::make_unique<SUNLinSolBand>(xB, model->ubw, model->lbw);
        setLinearSolverB(which);
        setBandJacFnB(which);
        break;

    case LinearSolver::LAPACKDense:
        throw AmiException("Solver currently not supported!");

    case LinearSolver::LAPACKBand:
        throw AmiException("Solver currently not supported!");

    case LinearSolver::diag:
        diagB(which);
        setDenseJacFnB(which);
        break;

        /* ITERATIVE SOLVERS */

    case LinearSolver::SPGMR:
        linearSolverB = std::make_unique<SUNLinSolSPGMR>(xB);
        setLinearSolverB(which);
        setJacTimesVecFnB(which);
        break;

    case LinearSolver::SPBCG:
        linearSolverB = std::make_unique<SUNLinSolSPBCGS>(xB);
        setLinearSolverB(which);
        setJacTimesVecFnB(which);
        break;

    case LinearSolver::SPTFQMR:
        linearSolverB = std::make_unique<SUNLinSolSPTFQMR>(xB);
        setLinearSolverB(which);
        setJacTimesVecFnB(which);
        break;

        /* SPARSE SOLVERS */

    case LinearSolver::KLU:
        linearSolverB = std::make_unique<SUNLinSolKLU>(
            xB, model->nnz, CSC_MAT,
            static_cast<SUNLinSolKLU::StateOrdering>(getStateOrdering()));
        setLinearSolverB(which);
        setSparseJacFnB(which);
        break;
#ifdef SUNDIALS_SUPERLUMT
    case LinearSolver::SuperLUMT:
        linearSolverB = std::make_unique<SUNLinSolSuperLUMT>(
            *xB, model->nnz, CSC_MAT,
            static_cast<SUNLinSolSuperLUMT::StateOrdering>(getStateOrdering()));
        setLinearSolverB(which);
        setSparseJacFnB(which);
        break;
#endif
    default:
        throw AmiException("Invalid choice of solver: %d",
                           static_cast<int>(linsol));
    }
}

void Solver::initializeNonLinearSolverB(const int which) const {
    switch (iter) {
    case NonlinearSolverIteration::newton:
        nonLinearSolverB =
            std::make_unique<SUNNonLinSolNewton>(xB.getNVector());
        break;
    case NonlinearSolverIteration::fixedpoint:
        nonLinearSolverB =
            std::make_unique<SUNNonLinSolFixedPoint>(xB.getNVector());
        break;
    default:
        throw AmiException("Invalid non-linear solver specified (%d).",
                           static_cast<int>(iter));
    }

    setNonLinearSolverB(which);
}

bool operator==(const Solver &a, const Solver &b) {
    if (typeid(a) != typeid(b))
        return false;

    return (a.interpType == b.interpType) && (a.lmm == b.lmm) &&
           (a.iter == b.iter) && (a.stldet == b.stldet) &&
           (a.ordering == b.ordering) &&
           (a.newton_maxsteps == b.newton_maxsteps) &&
           (a.newton_maxlinsteps == b.newton_maxlinsteps) &&
           (a.newton_damping_factor_mode == b.newton_damping_factor_mode) &&
           (a.newton_damping_factor_lower_bound == b.newton_damping_factor_lower_bound) &&
           (a.requires_preequilibration == b.requires_preequilibration) && (a.ism == b.ism) &&
           (a.linsol == b.linsol) && (a.atol == b.atol) && (a.rtol == b.rtol) &&
           (a.maxsteps == b.maxsteps) && (a.maxstepsB == b.maxstepsB) &&
           (a.quad_atol == b.quad_atol) && (a.quad_rtol == b.quad_rtol) &&
           (a.getAbsoluteToleranceSteadyState() ==
            b.getAbsoluteToleranceSteadyState()) &&
           (a.getRelativeToleranceSteadyState() ==
            b.getRelativeToleranceSteadyState()) &&
           (a.getAbsoluteToleranceSteadyStateSensi() ==
            b.getAbsoluteToleranceSteadyStateSensi()) &&
           (a.getRelativeToleranceSteadyStateSensi() ==
            b.getRelativeToleranceSteadyStateSensi()) &&
           (a.rtol_fsa == b.rtol_fsa ||
            (isNaN(a.rtol_fsa) && isNaN(b.rtol_fsa))) &&
           (a.atol_fsa == b.atol_fsa ||
            (isNaN(a.atol_fsa) && isNaN(b.atol_fsa))) &&
           (a.rtolB == b.rtolB || (isNaN(a.rtolB) && isNaN(b.rtolB))) &&
           (a.atolB == b.atolB || (isNaN(a.atolB) && isNaN(b.atolB))) &&
           (a.sensi == b.sensi) && (a.sensi_meth == b.sensi_meth);
}

void Solver::applyTolerances() const {
    if (!getInitDone())
        throw AmiException(("Solver instance was not yet set up, the "
                            "tolerances cannot be applied yet!"));

    setSStolerances(this->rtol, this->atol);
}

void Solver::applyTolerancesFSA() const {
    if (!getInitDone())
        throw AmiException(("Solver instance was not yet set up, the "
                            "tolerances cannot be applied yet!"));

    if (sensi < SensitivityOrder::first)
        return;

    if (nplist()) {
        std::vector<realtype> atols(nplist(), getAbsoluteToleranceFSA());
        setSensSStolerances(getRelativeToleranceFSA(), atols.data());
        setSensErrCon(true);
    }
}

void Solver::applyTolerancesASA(const int which) const {
    if (!getAdjInitDone())
        throw AmiException(("Adjoint solver instance was not yet set up, the "
                            "tolerances cannot be applied yet!"));

    if (sensi < SensitivityOrder::first)
        return;

    /* specify integration tolerances for backward problem */
    setSStolerancesB(which, getRelativeToleranceB(), getAbsoluteToleranceB());
}

void Solver::applyQuadTolerancesASA(const int which) const {
    if (!getAdjInitDone())
        throw AmiException(("Adjoint solver instance was not yet set up, the "
                            "tolerances cannot be applied yet!"));

    if (sensi < SensitivityOrder::first)
        return;

    realtype quad_rtol = isNaN(this->quad_rtol) ? rtol : this->quad_rtol;
    realtype quad_atol = isNaN(this->quad_atol) ? atol : this->quad_atol;

    /* Enable Quadrature Error Control */
    setQuadErrConB(which, !std::isinf(quad_atol) && !std::isinf(quad_rtol));

    quadSStolerancesB(which, quad_rtol, quad_atol);
}

void Solver::applySensitivityTolerances() const {
    if (sensi < SensitivityOrder::first)
        return;

    if (sensi_meth == SensitivityMethod::forward)
        applyTolerancesFSA();
    else if (sensi_meth == SensitivityMethod::adjoint && getAdjInitDone()) {
        for (int iMem = 0; iMem < (int)solverMemoryB.size(); ++iMem)
            applyTolerancesASA(iMem);
    }
}

SensitivityMethod Solver::getSensitivityMethod() const { return sensi_meth; }

void Solver::setSensitivityMethod(const SensitivityMethod sensi_meth) {
    if (sensi_meth != this->sensi_meth)
        resetMutableMemory(nx(), nplist(), nquad());
    this->sensi_meth = sensi_meth;
}

bool Solver::getNewtonSolverBackward() const { return newton_solver_backward; }

void Solver::setNewtonSolverBackward(bool newton_solver_backward) {
  this->newton_solver_backward = newton_solver_backward;
}

int Solver::getNewtonMaxSteps() const { return newton_maxsteps; }

void Solver::setNewtonMaxSteps(const int newton_maxsteps) {
    if (newton_maxsteps < 0)
        throw AmiException("newton_maxsteps must be a non-negative number");
    this->newton_maxsteps = newton_maxsteps;
}

bool Solver::getPreequilibration() const { return requires_preequilibration; }

void Solver::setPreequilibration(const bool require_preequilibration) {
    this->requires_preequilibration = require_preequilibration;
}

int Solver::getNewtonMaxLinearSteps() const { return newton_maxlinsteps; }

void Solver::setNewtonMaxLinearSteps(const int newton_maxlinsteps) {
    if (newton_maxlinsteps < 0)
        throw AmiException("newton_maxlinsteps must be a non-negative number");
    this->newton_maxlinsteps = newton_maxlinsteps;
}

NewtonDampingFactorMode Solver::getNewtonDampingFactorMode() const { return newton_damping_factor_mode; }

void Solver::setNewtonDampingFactorMode(NewtonDampingFactorMode dampingFactorMode) {
  this->newton_damping_factor_mode = dampingFactorMode;
}

double Solver::getNewtonDampingFactorLowerBound() const { return newton_damping_factor_lower_bound; }

void Solver::setNewtonDampingFactorLowerBound(double dampingFactorLowerBound) {
  this->newton_damping_factor_lower_bound = dampingFactorLowerBound;
}

SensitivityOrder Solver::getSensitivityOrder() const { return sensi; }

void Solver::setSensitivityOrder(const SensitivityOrder sensi) {
    if (this->sensi != sensi)
        resetMutableMemory(nx(), nplist(), nquad());
    this->sensi = sensi;

    if (getInitDone())
        applySensitivityTolerances();
}

double Solver::getRelativeTolerance() const {
    return static_cast<double>(rtol);
}

void Solver::setRelativeTolerance(const double rtol) {
    if (rtol < 0)
        throw AmiException("rtol must be a non-negative number");

    this->rtol = static_cast<realtype>(rtol);

    if (getInitDone()) {
        applyTolerances();
        applySensitivityTolerances();
    }
}

double Solver::getAbsoluteTolerance() const {
    return static_cast<double>(atol);
}

void Solver::setAbsoluteTolerance(double atol) {
    if (atol < 0)
        throw AmiException("atol must be a non-negative number");

    this->atol = static_cast<realtype>(atol);

    if (getInitDone()) {
        applyTolerances();
        applySensitivityTolerances();
    }
}

double Solver::getRelativeToleranceFSA() const {
    return static_cast<double>(isNaN(rtol_fsa) ? rtol : rtol_fsa);
}

void Solver::setRelativeToleranceFSA(const double rtol) {
    if (rtol < 0)
        throw AmiException("rtol must be a non-negative number");

    rtol_fsa = static_cast<realtype>(rtol);

    if (getInitDone()) {
        applySensitivityTolerances();
    }
}

double Solver::getAbsoluteToleranceFSA() const {
    return static_cast<double>(isNaN(atol_fsa) ? atol : atol_fsa);
}

void Solver::setAbsoluteToleranceFSA(const double atol) {
    if (atol < 0)
        throw AmiException("atol must be a non-negative number");

    atol_fsa = static_cast<realtype>(atol);

    if (getInitDone()) {
        applySensitivityTolerances();
    }
}

double Solver::getRelativeToleranceB() const {
    return static_cast<double>(isNaN(rtolB) ? rtol : rtolB);
}

void Solver::setRelativeToleranceB(const double rtol) {
    if (rtol < 0)
        throw AmiException("rtol must be a non-negative number");

    rtolB = static_cast<realtype>(rtol);

    if (getInitDone()) {
        applySensitivityTolerances();
    }
}

double Solver::getAbsoluteToleranceB() const {
    return static_cast<double>(isNaN(atolB) ? atol : atolB);
}

void Solver::setAbsoluteToleranceB(const double atol) {
    if (atol < 0)
        throw AmiException("atol must be a non-negative number");

    atolB = static_cast<realtype>(atol);

    if (getInitDone()) {
        applySensitivityTolerances();
    }
}

double Solver::getRelativeToleranceQuadratures() const {
    return static_cast<double>(quad_rtol);
}

void Solver::setRelativeToleranceQuadratures(const double rtol) {
    if (rtol < 0)
        throw AmiException("rtol must be a non-negative number");

    this->quad_rtol = static_cast<realtype>(rtol);

    if (sensi_meth != SensitivityMethod::adjoint)
        return;

    for (int iMem = 0; iMem < (int)solverMemoryB.size(); ++iMem)
        if (solverMemoryB.at(iMem))
            applyQuadTolerancesASA(iMem);
}

double Solver::getAbsoluteToleranceQuadratures() const {
    return static_cast<double>(quad_atol);
}

void Solver::setAbsoluteToleranceQuadratures(const double atol) {
    if (atol < 0)
        throw AmiException("atol must be a non-negative number");

    this->quad_atol = static_cast<realtype>(atol);

    if (sensi_meth != SensitivityMethod::adjoint)
        return;

    for (int iMem = 0; iMem < (int)solverMemoryB.size(); ++iMem)
        if (solverMemoryB.at(iMem))
            applyTolerancesASA(iMem);
}

double Solver::getRelativeToleranceSteadyState() const {
    return static_cast<double>(isNaN(ss_rtol) ? rtol : ss_rtol);
}

void Solver::setRelativeToleranceSteadyState(const double rtol) {
    if (rtol < 0)
        throw AmiException("rtol must be a non-negative number");

    this->ss_rtol = static_cast<realtype>(rtol);
}

double Solver::getAbsoluteToleranceSteadyState() const {
    return static_cast<double>(isNaN(ss_atol) ? atol : ss_atol);
}

void Solver::setAbsoluteToleranceSteadyState(const double atol) {
    if (atol < 0)
        throw AmiException("atol must be a non-negative number");

    this->ss_atol = static_cast<realtype>(atol);
}

double Solver::getRelativeToleranceSteadyStateSensi() const {
    return static_cast<double>(isNaN(ss_rtol_sensi) ? rtol : ss_rtol_sensi);
}

void Solver::setRelativeToleranceSteadyStateSensi(const double rtol) {
    if (rtol < 0)
        throw AmiException("rtol must be a non-negative number");

    this->ss_rtol_sensi = static_cast<realtype>(rtol);
}

double Solver::getAbsoluteToleranceSteadyStateSensi() const {
    return static_cast<double>(isNaN(ss_atol_sensi) ? atol : ss_atol_sensi);
}

void Solver::setAbsoluteToleranceSteadyStateSensi(const double atol) {
    if (atol < 0)
        throw AmiException("atol must be a non-negative number");

    this->ss_atol_sensi = static_cast<realtype>(atol);
}

long int Solver::getMaxSteps() const { return maxsteps; }

void Solver::setMaxSteps(const long int maxsteps) {
    if (maxsteps < 0)
        throw AmiException("maxsteps must be a non-negative number");

    this->maxsteps = maxsteps;
    if (getAdjInitDone())
        resetMutableMemory(nx(), nplist(), nquad());
}

long int Solver::getMaxStepsBackwardProblem() const { return maxstepsB; }

void Solver::setMaxStepsBackwardProblem(const long int maxsteps) {
    if (maxsteps < 0)
        throw AmiException("maxsteps must be a non-negative number");

    this->maxstepsB = maxsteps;
    for (int iMem = 0; iMem < (int)solverMemoryB.size(); ++iMem)
        if (solverMemoryB.at(iMem))
            setMaxNumStepsB(iMem, this->maxstepsB);
}

LinearMultistepMethod Solver::getLinearMultistepMethod() const { return lmm; }

void Solver::setLinearMultistepMethod(const LinearMultistepMethod lmm) {
    if (solverMemory)
        resetMutableMemory(nx(), nplist(), nquad());
    this->lmm = lmm;
}

NonlinearSolverIteration Solver::getNonlinearSolverIteration() const {
    return iter;
}

void Solver::setNonlinearSolverIteration(const NonlinearSolverIteration iter) {
    if (solverMemory)
        resetMutableMemory(nx(), nplist(), nquad());
    this->iter = iter;
}

InterpolationType Solver::getInterpolationType() const { return interpType; }

void Solver::setInterpolationType(const InterpolationType interpType) {
    if (!solverMemoryB.empty())
        resetMutableMemory(nx(), nplist(), nquad());
    this->interpType = interpType;
}

int Solver::getStateOrdering() const { return ordering; }

void Solver::setStateOrdering(int ordering) {
    this->ordering = ordering;
    if (solverMemory && linsol == LinearSolver::KLU) {
        auto klu = dynamic_cast<SUNLinSolKLU *>(linearSolver.get());
        klu->setOrdering(static_cast<SUNLinSolKLU::StateOrdering>(ordering));
        klu = dynamic_cast<SUNLinSolKLU *>(linearSolverB.get());
        klu->setOrdering(static_cast<SUNLinSolKLU::StateOrdering>(ordering));
    }
#ifdef SUNDIALS_SUPERLUMT
    if (solverMemory && linsol == LinearSolver::SuperLUMT) {
        auto klu = dynamic_cast<SUNLinSolSuperLUMT *>(linearSolver.get());
        klu->setOrdering(
            static_cast<SUNLinSolSuperLUMT::StateOrdering>(ordering));
        klu = dynamic_cast<SUNLinSolSuperLUMT *>(linearSolverB.get());
        klu->setOrdering(
            static_cast<SUNLinSolSuperLUMT::StateOrdering>(ordering));
    }
#endif
}

int Solver::getStabilityLimitFlag() const { return stldet; }

void Solver::setStabilityLimitFlag(const int stldet) {
    if (stldet != TRUE && stldet != FALSE)
        throw AmiException("Invalid stldet flag, valid values are %i or %i",
                           TRUE, FALSE);

    this->stldet = stldet;
    if (solverMemory) {
        setStabLimDet(stldet);
        for (int iMem = 0; iMem < (int)solverMemoryB.size(); ++iMem)
            if (solverMemoryB.at(iMem))
                setStabLimDetB(iMem, stldet);
    }
}

LinearSolver Solver::getLinearSolver() const { return linsol; }

void Solver::setLinearSolver(LinearSolver linsol) {
    if (solverMemory)
        resetMutableMemory(nx(), nplist(), nquad());
    this->linsol = linsol;
    /*if(solverMemory)
     initializeLinearSolver(getModel());*/
}

InternalSensitivityMethod Solver::getInternalSensitivityMethod() const {
    return ism;
}

void Solver::setInternalSensitivityMethod(const InternalSensitivityMethod ism) {
    if (solverMemory)
        resetMutableMemory(nx(), nplist(), nquad());
    this->ism = ism;
}

void Solver::initializeNonLinearSolverSens(const Model *model) const {
    switch (iter) {
    case NonlinearSolverIteration::newton:
        switch (ism) {
        case InternalSensitivityMethod::staggered:
        case InternalSensitivityMethod::simultaneous:
            nonLinearSolverSens = std::make_unique<SUNNonLinSolNewton>(
                1 + model->nplist(), x.getNVector());
            break;
        case InternalSensitivityMethod::staggered1:
            nonLinearSolverSens =
                std::make_unique<SUNNonLinSolNewton>(x.getNVector());
            break;
        default:
            throw AmiException(
                "Unsupported internal sensitivity method selected: %d", ism);
        }
        break;
    case NonlinearSolverIteration::fixedpoint:
        switch (ism) {
        case InternalSensitivityMethod::staggered:
        case InternalSensitivityMethod::simultaneous:
            nonLinearSolverSens = std::make_unique<SUNNonLinSolFixedPoint>(
                1 + model->nplist(), x.getNVector());
            break;
        case InternalSensitivityMethod::staggered1:
            nonLinearSolverSens =
                std::make_unique<SUNNonLinSolFixedPoint>(x.getNVector());
            break;
        default:
            throw AmiException(
                "Unsupported internal sensitivity method selected: %d", ism);
        }
        break;
    default:
        throw AmiException("Invalid non-linear solver specified (%d).",
                           static_cast<int>(iter));
    }

    setNonLinearSolverSens();
}

int Solver::nplist() const { return sx.getLength(); }

int Solver::nx() const { return x.getLength(); }

int Solver::nquad() const { return xQB.getLength(); }

bool Solver::getInitDone() const { return initialized; };

bool Solver::getSensInitDone() const { return sensInitialized; }

bool Solver::getAdjInitDone() const { return adjInitialized; }

bool Solver::getInitDoneB(const int which) const {
    return static_cast<int>(initializedB.size()) > which &&
           initializedB.at(which);
}

bool Solver::getQuadInitDoneB(const int which) const {
    return static_cast<int>(initializedQB.size()) > which &&
           initializedQB.at(which);
}

void Solver::setInitDone() const { initialized = true; };

void Solver::setSensInitDone() const { sensInitialized = true; }

void Solver::setAdjInitDone() const { adjInitialized = true; }

void Solver::setInitDoneB(const int which) const {
    if (which >= static_cast<int>(initializedB.size()))
        initializedB.resize(which + 1, false);
    initializedB.at(which) = true;
}

void Solver::setQuadInitDoneB(const int which) const {
    if (which >= static_cast<int>(initializedQB.size()))
        initializedQB.resize(which + 1, false);
    initializedQB.at(which) = true;
}

realtype Solver::getCpuTime() const {
    return cpu_time;
}

realtype Solver::getCpuTimeB() const {
    return cpu_timeB;
}

void Solver::resetMutableMemory(const int nx, const int nplist,
                                const int nquad) const {
    solverMemory = nullptr;
    initialized = false;
    adjInitialized = false;
    sensInitialized = false;
    solverWasCalledF = false;
    solverWasCalledB = false;

    x = AmiVector(nx);
    dx = AmiVector(nx);
    sx = AmiVectorArray(nx, nplist);
    sdx = AmiVectorArray(nx, nplist);

    xB = AmiVector(nx);
    dxB = AmiVector(nx);
    xQB = AmiVector(nquad);

    solverMemoryB.clear();
    initializedB.clear();
    initializedQB.clear();
}

void Solver::writeSolution(realtype *t, AmiVector &x, AmiVector &dx,
                           AmiVectorArray &sx) const {
    *t = gett();
    x.copy(getState(*t));
    dx.copy(getDerivativeState(*t));
    if (sensInitialized) {
        sx.copy(getStateSensitivity(*t));
    }
}

void Solver::writeSolutionB(realtype *t, AmiVector &xB, AmiVector &dxB,
                            AmiVector &xQB, const int which) const {
    *t = gett();
    xB.copy(getAdjointState(which, *t));
    dxB.copy(getAdjointDerivativeState(which, *t));
    xQB.copy(getAdjointQuadrature(which, *t));
}

const AmiVector &Solver::getState(const realtype t) const {
    if (t == this->t)
        return x;

    if (solverWasCalledF)
        getDky(t, 0);

    return dky;
}

const AmiVector &Solver::getDerivativeState(const realtype t) const {
    if (t == this->t)
        return dx;

    if (solverWasCalledF)
        getDky(t, 1);

    return dky;
}

const AmiVectorArray &Solver::getStateSensitivity(const realtype t) const {
    if (sensInitialized) {
        if (solverWasCalledF) {
            if (t == this->t) {
                getSens();
            } else {
                getSensDky(t, 0);
            }
        }
    } else {
        sx.reset();
    }
    return sx;
}

const AmiVector &Solver::getAdjointState(const int which,
                                         const realtype t) const {
    if (adjInitialized) {
        if (solverWasCalledB) {
            if (t == this->t) {
                getB(which);
                return xB;
            }
            getDkyB(t, 0, which);
        }
    } else {
        dky.reset();
    }
    return dky;
}

const AmiVector &Solver::getAdjointDerivativeState(const int which,
                                                   const realtype t) const {
    if (adjInitialized) {
        if (solverWasCalledB) {
            if (t == this->t) {
                getB(which);
                return dxB;
            }
            getDkyB(t, 1, which);
        }
    } else {
        dky.reset();
    }
    return dky;
}

const AmiVector &Solver::getAdjointQuadrature(const int which,
                                              const realtype t) const {
    if (adjInitialized) {
        if (solverWasCalledB) {
            if (t == this->t) {
                getQuadB(which);
                return xQB;
            }
            getQuadDkyB(t, 0, which);
        }
    } else {
        xQB.reset();
    }
    return xQB;
}

realtype Solver::gett() const { return t; }

void wrapErrHandlerFn(int error_code, const char *module,
                      const char *function, char *msg, void * eh_data) {
#define BUF_SIZE 250
    char buffer[BUF_SIZE];
    char buffid[BUF_SIZE];
    snprintf(buffer, BUF_SIZE, "AMICI ERROR: in module %s in function %s : %s ", module,
            function, msg);
    switch (error_code) {
    case 99:
        snprintf(buffid, BUF_SIZE, "AMICI:%s:%s:WARNING", module, function);
        break;

    case -1:
        snprintf(buffid, BUF_SIZE, "AMICI:%s:%s:TOO_MUCH_WORK", module, function);
        break;

    case -2:
        snprintf(buffid, BUF_SIZE, "AMICI:%s:%s:TOO_MUCH_ACC", module, function);
        break;

    case -3:
        snprintf(buffid, BUF_SIZE, "AMICI:%s:%s:ERR_FAILURE", module, function);
        break;

    case -4:
        snprintf(buffid, BUF_SIZE, "AMICI:%s:%s:CONV_FAILURE", module, function);
        break;

    default:
        snprintf(buffid, BUF_SIZE, "AMICI:%s:%s:OTHER", module, function);
        break;
    }


    if(!eh_data) {
        throw std::runtime_error("eh_data unset");
    }
    auto solver = static_cast<Solver const*>(eh_data);
    solver->app->warning(buffid, buffer);
}

} // namespace amici
