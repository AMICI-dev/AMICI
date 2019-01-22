#include <amici/sundials_linsol_wrapper.h>
#include <amici/exception.h>

#include <utility>
#include <new> // bad_alloc

namespace amici {

SUNLinSolWrapper::SUNLinSolWrapper(SUNLinearSolver linsol_)
{
    linsol = linsol_;
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

SUNLinSolWrapper &SUNLinSolWrapper::operator=(SUNLinSolWrapper &&other) noexcept
{
    return *this = SUNLinSolWrapper(other.linsol);
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
    return SUNLinSolInitialize(linsol);
}

int SUNLinSolWrapper::setup(SUNMatrix A)
{
    return SUNLinSolSetup(linsol, A);
}

int SUNLinSolWrapper::setup(SUNMatrixWrapper A)
{
    return SUNLinSolSetup(linsol, A.get());
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





} // namespace amici
