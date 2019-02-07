#include "amici/solver.h"
#include "amici/exception.h"
#include "amici/forwardproblem.h"
#include "amici/backwardproblem.h"
#include "amici/model.h"
#include "amici/rdata.h"

#include <cmath>
#include <cstdio>
#include <cstring>

#include <sundials/sundials_spgmr.h>
#include <cvodes/cvodes_spils.h>

namespace amici {

extern msgIdAndTxtFp warnMsgIdAndTxt;


Solver::Solver(const Solver &other) : Solver()
{
    sensi = other.sensi;
    atol = other.atol;
    rtol = other.rtol;
    atol_fsa = other.atol_fsa;
    rtol_fsa = other.rtol_fsa;
    atolB = other.atolB;
    rtolB = other.rtolB;
    quad_atol = other.quad_atol;
    quad_rtol = other.quad_rtol;
    ss_atol = other.ss_atol;
    ss_rtol = other.ss_rtol;
    ss_atol_sensi = other.ss_atol_sensi;
    ss_rtol_sensi = other.ss_rtol_sensi;
    maxsteps = other.maxsteps;
    maxstepsB = other.maxstepsB;
    newton_maxsteps = other.newton_maxsteps;
    newton_maxlinsteps = other.newton_maxlinsteps;
    newton_preeq = other.newton_preeq;
    ism = other.ism;
    sensi_meth = other.sensi_meth;
    linsol = other.linsol;
    interpType = other.interpType;
    lmm = other.lmm;
    iter = other.iter;
    stldet = other.stldet;
    ordering = other.ordering;
}

void Solver::setup(AmiVector *x, AmiVector *dx, AmiVectorArray *sx, AmiVectorArray *sdx, Model *model) {
    bool computeSensitivities = sensi >= SensitivityOrder::first
            && model->nx_solver > 0;
    model->initialize(x, dx, sx, sdx, computeSensitivities);

    /* Create solver memory object */
    allocateSolver();
    if (!solverMemory)
        throw AmiException("Failed to allocated solver memory!");

    /* Initialize CVodes/IDAs solver*/
    init(x, dx, model->t0());

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

    if (computeSensitivities) {
        auto plist = model->getParameterList();

        if (sensi_meth == SensitivityMethod::forward && !plist.empty()) {
            /* Set sensitivity analysis optional inputs */
            auto par = model->getUnscaledParameters();

            /* Activate sensitivity calculations */
            sensInit1(sx, sdx, plist.size());
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
    if(model->nt()>1)
        calcIC(model->t(1), x, dx);
}

void Solver::setupB(BackwardProblem *bwd, Model *model) {
    if (!solverMemory)
        throw AmiException("Solver for the forward problem must be setup first");

    /* write initial conditions */
    std::vector<realtype> dJydx = bwd->getdJydx();
    AmiVector *xB = bwd->getxBptr();
    xB->reset();
    for (int ix = 0; ix < model->nxtrue_solver; ++ix)
        for (int iJ = 0; iJ < model->nJ; ++iJ)
            xB->at(ix + iJ * model->nxtrue_solver) +=
                dJydx.at(iJ + ( ix + (model->nt() - 1)  * model->nx_solver ) * model->nJ);
    bwd->getdxBptr()->reset();
    bwd->getxQBptr()->reset();

    /* allocate memory for the backward problem */
    allocateSolverB(bwd->getwhichptr());

    /* initialise states */
    binit(bwd->getwhich(), bwd->getxBptr(), bwd->getdxBptr(), bwd->gett());

    /* Attach user data */
    setUserDataB(bwd->getwhich(), model);

    /* Number of maximal internal steps */
    setMaxNumStepsB(bwd->getwhich(), (maxstepsB == 0) ? maxsteps * 100 : maxstepsB);

    initializeLinearSolverB(model, bwd->getwhich());

    /* Initialise quadrature calculation */
    qbinit(bwd->getwhich(), bwd->getxQBptr());

    applyTolerancesASA(bwd->getwhich());
    applyQuadTolerancesASA(bwd->getwhich());

    setStabLimDetB(bwd->getwhich(), stldet);
}

void Solver::wrapErrHandlerFn(int error_code, const char *module,
                              const char *function, char *msg, void * /*eh_data*/) {
    char buffer[250];
    char buffid[250];
    sprintf(buffer, "AMICI ERROR: in module %s in function %s : %s ", module,
            function, msg);
    switch (error_code) {
    case 99:
        sprintf(buffid, "AMICI:mex:%s:%s:WARNING", module, function);
        break;

    case -1:
        sprintf(buffid, "AMICI:mex:%s:%s:TOO_MUCH_WORK", module, function);
        break;

    case -2:
        sprintf(buffid, "AMICI:mex:%s:%s:TOO_MUCH_ACC", module, function);
        break;

    case -3:
        sprintf(buffid, "AMICI:mex:%s:%s:ERR_FAILURE", module, function);
        break;

    case -4:
        sprintf(buffid, "AMICI:mex:%s:%s:CONV_FAILURE", module, function);
        break;

    default:
        sprintf(buffid, "AMICI:mex:%s:%s:OTHER", module, function);
        break;
    }

    warnMsgIdAndTxt(buffid, buffer);
}

void Solver::getDiagnosis(const int it, ReturnData *rdata) const {
    long int number;

    if(solverWasCalled && solverMemory) {
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

void Solver::getDiagnosisB(const int it, ReturnData *rdata, int which) const {
    long int number;

    if(solverWasCalled && solverMemoryB.at(which)) {
        getNumSteps(solverMemoryB.at(which).get(), &number);
        rdata->numstepsB[it] = (double)number;

        getNumRhsEvals(solverMemoryB.at(which).get(), &number);
        rdata->numrhsevalsB[it] = (double)number;

        getNumErrTestFails(solverMemoryB.at(which).get(), &number);
        rdata->numerrtestfailsB[it] = (double)number;

        getNumNonlinSolvConvFails(solverMemoryB.at(which).get(), &number);
        rdata->numnonlinsolvconvfailsB[it] = (double)number;
    }
}

/**
 * initializeLinearSolver sets the linear solver for the forward problem
 *
 * @param model pointer to the model object
 */
void Solver::initializeLinearSolver(const Model *model) {
    /* Attach linear solver module */

    switch (linsol) {

            /* DIRECT SOLVERS */

        case LinearSolver::dense:
            dense(model->nx_solver);
            setDenseJacFn();
            break;

        case LinearSolver::band:
            band(model->nx_solver, model->ubw, model->lbw);
            setBandJacFn();
            break;

        case LinearSolver::LAPACKDense:
            throw AmiException("Solver currently not supported!");
            /* status = CVLapackDense(ami_mem, nx_solver);
             if (status != AMICI_SUCCESS) return;

             status = SetDenseJacFn(ami_mem);
             if (status != AMICI_SUCCESS) return;
             */

        case LinearSolver::LAPACKBand:
            throw AmiException("Solver currently not supported!");
            /* status = CVLapackBand(ami_mem, nx_solver);
             if (status != AMICI_SUCCESS) return;

             status = SetBandJacFn(ami_mem);
             if (status != AMICI_SUCCESS) return;
             */

        case LinearSolver::diag:
            diag();
            break;


            /* ITERATIVE SOLVERS */

        case LinearSolver::SPGMR:
            spgmr(PREC_NONE, CVSPILS_MAXL);
            setJacTimesVecFn();
            break;

        case LinearSolver::SPBCG:
            spbcg(PREC_NONE, CVSPILS_MAXL);
            setJacTimesVecFn();
            break;

        case LinearSolver::SPTFQMR:
            sptfqmr(PREC_NONE, CVSPILS_MAXL);
            setJacTimesVecFn();
            break;

            /* SPARSE SOLVERS */

        case LinearSolver::KLU:
            klu(model->nx_solver, model->nnz, CSC_MAT);
            setSparseJacFn();
            kluSetOrdering((int) getStateOrdering());
            break;

        default:
            throw AmiException("Invalid choice of solver!");

    }
}

void Solver::initializeLinearSolverB(const Model *model, const int which) {
    switch (linsol) {

            /* DIRECT SOLVERS */

        case LinearSolver::dense:
            denseB(which, model->nx_solver);
            setDenseJacFnB(which);
            break;

        case LinearSolver::band:
            bandB(which, model->nx_solver, model->ubw, model->lbw);
            setBandJacFnB(which);
            break;

        case LinearSolver::LAPACKDense:

            /* #if SUNDIALS_BLAS_LAPACK
             status = CVLapackDenseB(ami_mem, bwd->getwhich(), nx_solver);
             if (status != AMICI_SUCCESS) return;

             status = SetDenseJacFnB(ami_mem, bwd->getwhich());
             if (status != AMICI_SUCCESS) return;
             #else*/
            throw AmiException("Solver currently not supported!");
            /* #endif*/

        case LinearSolver::LAPACKBand:

            /* #if SUNDIALS_BLAS_LAPACK
             status = CVLapackBandB(ami_mem, bwd->getwhich(), nx_solver, ubw, lbw);
             if (status != AMICI_SUCCESS) return;

             status = SetBandJacFnB(ami_mem, bwd->getwhich());
             if (status != AMICI_SUCCESS) return;
             #else*/
            throw AmiException("Solver currently not supported!");
            /* #endif*/

        case LinearSolver::diag:
            diagB(which);
            setDenseJacFnB(which);
            break;

            /* ITERATIVE SOLVERS */

        case LinearSolver::SPGMR:
            spgmrB(which, PREC_NONE, CVSPILS_MAXL);
            setJacTimesVecFnB(which);
            break;

        case LinearSolver::SPBCG:
            spbcgB(which, PREC_NONE, CVSPILS_MAXL);
            setJacTimesVecFnB(which);
            break;

        case LinearSolver::SPTFQMR:
            sptfqmrB(which, PREC_NONE, CVSPILS_MAXL);
            setJacTimesVecFnB(which);
            break;

            /* SPARSE SOLVERS */

        case LinearSolver::KLU:
            kluB(which, model->nx_solver, model->nnz, CSC_MAT);
            setSparseJacFnB(which);
            kluSetOrderingB(which, (int) getStateOrdering());
            break;

        default:
            throw AmiException("Invalid local Solver!");
    }
}

bool operator ==(const Solver &a, const Solver &b)
{
    if (typeid(a) != typeid(b))
            return false;

    return (a.interpType == b.interpType)
            && (a.lmm == b.lmm)
            && (a.iter == b.iter)
            && (a.stldet == b.stldet)
            && (a.ordering == b.ordering)
            && (a.newton_maxsteps == b.newton_maxsteps)
            && (a.newton_maxlinsteps == b.newton_maxlinsteps)
            && (a.newton_preeq == b.newton_preeq)
            && (a.ism == b.ism)
            && (a.linsol == b.linsol)
            && (a.atol == b.atol)
            && (a.rtol == b.rtol)
            && (a.maxsteps == b.maxsteps)
            && (a.maxstepsB == b.maxstepsB)
            && (a.quad_atol == b.quad_atol)
            && (a.quad_rtol == b.quad_rtol)
            && (a.getAbsoluteToleranceSteadyState() ==
                b.getAbsoluteToleranceSteadyState())
            && (a.getRelativeToleranceSteadyState() ==
                b.getRelativeToleranceSteadyState())
            && (a.getAbsoluteToleranceSteadyStateSensi() ==
                b.getAbsoluteToleranceSteadyStateSensi())
            && (a.getRelativeToleranceSteadyStateSensi() ==
                b.getRelativeToleranceSteadyStateSensi())
            && (a.rtol_fsa == b.rtol_fsa || (isNaN(a.rtol_fsa) && isNaN(b.rtol_fsa)))
            && (a.atol_fsa == b.atol_fsa || (isNaN(a.atol_fsa) && isNaN(b.atol_fsa)))
            && (a.rtolB == b.rtolB || (isNaN(a.rtolB) && isNaN(b.rtolB)))
            && (a.atolB == b.atolB || (isNaN(a.atolB) && isNaN(b.atolB)))
            && (a.sensi == b.sensi)
            && (a.sensi_meth == b.sensi_meth);
}

void Solver::applyTolerances() {
    if (!getMallocDone())
        throw AmiException(("Solver instance was not yet set up, the tolerances cannot be applied yet!"));

    setSStolerances(RCONST(this->rtol), RCONST(this->atol));
}

void Solver::applyTolerancesFSA() {
    if (!getMallocDone())
        throw AmiException(("Solver instance was not yet set up, the tolerances cannot be applied yet!"));

    if (sensi < SensitivityOrder::first)
        return;

    if(nplist()) {
        std::vector<realtype> atols(nplist(), getAbsoluteToleranceFSA());
        setSensSStolerances(getRelativeToleranceFSA(), atols.data());
        setSensErrCon(true);
    }
}

void Solver::applyTolerancesASA(int which) {
    if (!getAdjMallocDone())
        throw AmiException(("Adjoint solver instance was not yet set up, the tolerances cannot be applied yet!"));

    if (sensi < SensitivityOrder::first)
        return;

    /* specify integration tolerances for backward problem */
    setSStolerancesB(which, RCONST(getRelativeToleranceB()), RCONST(getAbsoluteToleranceB()));
}

void Solver::applyQuadTolerancesASA(int which) {
    if (!getAdjMallocDone())
        throw AmiException(("Adjoint solver instance was not yet set up, the tolerances cannot be applied yet!"));

    if (sensi < SensitivityOrder::first)
        return;

    double quad_rtol = isNaN(this->quad_rtol) ? rtol : this->quad_rtol;
    double quad_atol = isNaN(this->quad_atol) ? atol : this->quad_atol;

    /* Enable Quadrature Error Control */
    setQuadErrConB(which,
                   !std::isinf(quad_atol) && !std::isinf(quad_rtol));

    quadSStolerancesB(which,
                      RCONST(quad_rtol),
                      RCONST(quad_atol));
}

void Solver::applySensitivityTolerances() {
    if (sensi < SensitivityOrder::first)
        return;

    if (sensi_meth == SensitivityMethod::forward)
        applyTolerancesFSA();
    else if (sensi_meth == SensitivityMethod::adjoint && getAdjMallocDone()) {
        for (int iMem = 0; iMem < (int) solverMemoryB.size(); ++iMem)
            applyTolerancesASA(iMem);
    }
}

SensitivityMethod Solver::getSensitivityMethod() const{
    return sensi_meth;
}

void Solver::setSensitivityMethod(SensitivityMethod sensi_meth) {
    this->sensi_meth = sensi_meth;
}

int Solver::getNewtonMaxSteps() const {
    return newton_maxsteps;
}

void Solver::setNewtonMaxSteps(int newton_maxsteps) {
    if(newton_maxsteps < 0)
        throw AmiException("newton_maxsteps must be a non-negative number");
    this->newton_maxsteps = newton_maxsteps;
}

bool Solver::getNewtonPreequilibration() const {
    return newton_preeq;
}

void Solver::setNewtonPreequilibration(bool newton_preeq) {
    this->newton_preeq = newton_preeq;
}

int Solver::getNewtonMaxLinearSteps() const {
    return newton_maxlinsteps;
}

void Solver::setNewtonMaxLinearSteps(int newton_maxlinsteps) {
    if(newton_maxlinsteps < 0)
        throw AmiException("newton_maxlinsteps must be a non-negative number");
    this->newton_maxlinsteps = newton_maxlinsteps;
}

SensitivityOrder Solver::getSensitivityOrder() const {
    return sensi;
}

void Solver::setSensitivityOrder(SensitivityOrder sensi) {
    this->sensi = sensi;

    if(getMallocDone())
        applySensitivityTolerances();
}

double Solver::getRelativeTolerance() const {
    return rtol;
}

void Solver::setRelativeTolerance(double rtol) {
    if(rtol < 0)
        throw AmiException("rtol must be a non-negative number");

    this->rtol = rtol;

    if(getMallocDone()) {
        applyTolerances();
        applySensitivityTolerances();
    }
}

double Solver::getAbsoluteTolerance() const {
    return atol;
}

void Solver::setAbsoluteTolerance(double atol) {
    if(atol < 0)
        throw AmiException("atol must be a non-negative number");

    this->atol = atol;

    if(getMallocDone()) {
        applyTolerances();
        applySensitivityTolerances();
    }
}

double Solver::getRelativeToleranceFSA() const {
    return isNaN(rtol_fsa) ? rtol : rtol_fsa;
}

void Solver::setRelativeToleranceFSA(double rtol) {
    if(rtol < 0)
        throw AmiException("rtol must be a non-negative number");

    rtol_fsa = rtol;

    if(getMallocDone()) {
        applySensitivityTolerances();
    }
}

double Solver::getAbsoluteToleranceFSA() const {
    return isNaN(atol_fsa) ? atol : atol_fsa;
}

void Solver::setAbsoluteToleranceFSA(double atol) {
    if(atol < 0)
        throw AmiException("atol must be a non-negative number");

    atol_fsa = atol;

    if(getMallocDone()) {
        applySensitivityTolerances();
    }
}

double Solver::getRelativeToleranceB() const
{
    return isNaN(rtolB) ? rtol : rtolB;
}

void Solver::setRelativeToleranceB(double rtol)
{
    if(rtol < 0)
        throw AmiException("rtol must be a non-negative number");

    rtolB = rtol;

    if(getMallocDone()) {
        applySensitivityTolerances();
    }
}

double Solver::getAbsoluteToleranceB() const
{
    return isNaN(atolB) ? atol : atolB;
}

void Solver::setAbsoluteToleranceB(double atol)
{
    if(atol < 0)
        throw AmiException("atol must be a non-negative number");

    atolB = atol;

    if(getMallocDone()) {
        applySensitivityTolerances();
    }
}

double Solver::getRelativeToleranceQuadratures() const {
    return quad_rtol;
}

void Solver::setRelativeToleranceQuadratures(double rtol) {
    if(rtol < 0)
        throw AmiException("rtol must be a non-negative number");

    this->quad_rtol = rtol;

    if (sensi_meth != SensitivityMethod::adjoint)
        return;

    for (int iMem = 0; iMem < (int) solverMemoryB.size(); ++iMem)
        if(solverMemoryB.at(iMem))
            applyQuadTolerancesASA(iMem);
}

double Solver::getAbsoluteToleranceQuadratures() const {
    return quad_atol;
}

void Solver::setAbsoluteToleranceQuadratures(double atol) {
    if (atol < 0)
        throw AmiException("atol must be a non-negative number");

    this->quad_atol = atol;

    if (sensi_meth != SensitivityMethod::adjoint)
        return;

    for (int iMem = 0; iMem < (int) solverMemoryB.size(); ++iMem)
        if(solverMemoryB.at(iMem))
            applyTolerancesASA(iMem);
}

double Solver::getRelativeToleranceSteadyState() const {
    return isNaN(ss_rtol) ? rtol : ss_rtol;
}

void Solver::setRelativeToleranceSteadyState(double rtol) {
    if(rtol < 0)
        throw AmiException("rtol must be a non-negative number");

    this->ss_rtol = rtol;
}

double Solver::getAbsoluteToleranceSteadyState() const {
    return isNaN(ss_atol) ? atol : ss_atol;
}

void Solver::setAbsoluteToleranceSteadyState(double atol) {
    if (atol < 0)
        throw AmiException("atol must be a non-negative number");

    this->ss_atol = atol;
}

double Solver::getRelativeToleranceSteadyStateSensi() const {
    return isNaN(ss_rtol_sensi) ? rtol : ss_rtol_sensi;
}

void Solver::setRelativeToleranceSteadyStateSensi(double rtol) {
    if(rtol < 0)
        throw AmiException("rtol must be a non-negative number");

    this->ss_rtol_sensi = rtol;
}

double Solver::getAbsoluteToleranceSteadyStateSensi() const {
    return isNaN(ss_atol_sensi) ? atol : ss_atol_sensi;
}

void Solver::setAbsoluteToleranceSteadyStateSensi(double atol) {
    if (atol < 0)
        throw AmiException("atol must be a non-negative number");

    this->ss_atol_sensi = atol;
}

int Solver::getMaxSteps() const {
    return maxsteps;
}

void Solver::setMaxSteps(int maxsteps) {
    if (maxsteps < 0)
        throw AmiException("maxsteps must be a non-negative number");

    this->maxsteps = maxsteps;
    if(solverMemory)
        setMaxNumSteps(this->maxsteps);
}

int Solver::getMaxStepsBackwardProblem() const {
    return maxstepsB;
}

void Solver::setMaxStepsBackwardProblem(int maxsteps) {
    if (maxsteps < 0)
        throw AmiException("maxsteps must be a non-negative number");

    this->maxstepsB = maxsteps;
    for (int iMem = 0; iMem < (int) solverMemoryB.size(); ++iMem)
        if(solverMemoryB.at(iMem))
            setMaxNumStepsB(iMem, this->maxstepsB);
}

LinearMultistepMethod Solver::getLinearMultistepMethod() const {
    return lmm;
}

void Solver::setLinearMultistepMethod(LinearMultistepMethod lmm) {
    if(solverMemory)
        throw AmiException("Solver instance was already set up, the linear system multistep method can no longer be changed!");
    this->lmm = lmm;
}

NonlinearSolverIteration Solver::getNonlinearSolverIteration() const {
    return iter;
}

void Solver::setNonlinearSolverIteration(NonlinearSolverIteration iter) {
    if(solverMemory)
        throw AmiException("Solver instance was already set up, the nonlinear system solution method can no longer be changed!");
    this->iter = iter;
}

InterpolationType Solver::getInterpolationType() const {
    return interpType;
}

void Solver::setInterpolationType(InterpolationType interpType) {
    if(!solverMemoryB.empty())
        throw AmiException("Adjoint solver object was already set up, the interpolation type can no longer be changed!");
    this->interpType = interpType;
}

StateOrdering Solver::getStateOrdering() const {
    return ordering;
}

void Solver::setStateOrdering(StateOrdering ordering) {
    this->ordering = ordering;
    if (solverMemory && linsol == LinearSolver::KLU) {
        kluSetOrdering((int)ordering);
        for (int iMem = 0; iMem < (int) solverMemoryB.size(); ++iMem)
            if(solverMemoryB.at(iMem))
                kluSetOrderingB(iMem, (int) ordering);
    }
}

int Solver::getStabilityLimitFlag() const {
    return stldet;
}

void Solver::setStabilityLimitFlag(int stldet) {
    if (stldet != TRUE && stldet != FALSE)
        throw AmiException("Invalid stldet flag, valid values are %i or %i",TRUE,FALSE);

    this->stldet = stldet;
    if (solverMemory) {
        setStabLimDet(stldet);
        for (int iMem = 0; iMem < (int) solverMemoryB.size(); ++iMem)
            if(solverMemoryB.at(iMem))
                setStabLimDetB(iMem,stldet);
    }
}

LinearSolver Solver::getLinearSolver() const {
    return linsol;
}

void Solver::setLinearSolver(LinearSolver linsol) {
    if(solverMemory)
        throw AmiException("Solver object was already set up, the linear solver can no longer be changed!");
    this->linsol = linsol;
    /*if(solverMemory)
     initializeLinearSolver(getModel());*/
}

InternalSensitivityMethod Solver::getInternalSensitivityMethod() const {
    return ism;
}

void Solver::setInternalSensitivityMethod(InternalSensitivityMethod ism) {
    if(solverMemory)
        throw AmiException("Solver object was already set up, the sensitivity method can no longer be changed!");
    this->ism = ism;
}


} // namespace amici
