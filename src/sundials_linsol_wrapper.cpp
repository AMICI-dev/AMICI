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
    return *this = SUNLinSolWrapper(linsol);
}

SUNLinearSolver SUNLinSolWrapper::get() const
{
    return linsol;
}

SUNLinearSolver_Type SUNLinSolWrapper::getType() const
{
    return SUNLinSolGetType(linsol);
}


} // namespace amici
