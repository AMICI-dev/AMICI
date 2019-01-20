#include <amici/sundials_matrix_wrapper.h>
#include <amici/exception.h>

#include <utility>
#include <new> // bad_alloc

namespace amici {

SUNMatrixWrapper::SUNMatrixWrapper(int M, int N, int NNZ, int sparsetype)
    : matrix(SUNSparseMatrix(M, N, NNZ, sparsetype))
{
    if(NNZ && !matrix)
        throw std::bad_alloc();
}
    
SUNMatrixWrapper::SUNMatrixWrapper(int M, int N)
: matrix(SUNDenseMatrix(M, N))
{
    if(M*N && !matrix)
        throw std::bad_alloc();
}

SUNMatrixWrapper::SUNMatrixWrapper(::SUNMatrix mat)
    : matrix(mat)
{}

SUNMatrixWrapper::~SUNMatrixWrapper()
{
    if(matrix)
        SUNMatDestroy(matrix);
}

SUNMatrixWrapper::SUNMatrixWrapper(const SUNMatrixWrapper &other)
{
    if(!other.matrix)
        return;

    matrix = SUNMatClone(other.matrix);
    if(!matrix)
        throw std::bad_alloc();

    SUNMatCopy(other.matrix, matrix);
}

SUNMatrixWrapper::SUNMatrixWrapper(SUNMatrixWrapper &&other) noexcept
{
    std::swap(matrix, other.matrix);

}

SUNMatrixWrapper &SUNMatrixWrapper::operator=(const SUNMatrixWrapper &other)
{
    return *this = SUNMatrixWrapper(other);
}

SUNMatrixWrapper &SUNMatrixWrapper::operator=(SUNMatrixWrapper &&other) noexcept
{
    std::swap(matrix, other.matrix);
    return *this;
}

realtype *SUNMatrixWrapper::data() {
    if(!matrix)
        return nullptr;
    
    switch(SUNMatGetID(matrix)) {
        case SUNMATRIX_DENSE:
            return SM_DATA_D(matrix);
        case SUNMATRIX_SPARSE:
            return SM_DATA_S(matrix);
        case SUNMATRIX_BAND:
            return SM_DATA_B(matrix);
        case SUNMATRIX_CUSTOM:
            throw AmiException("Amici currently does not support custom matrix "
                               "types.");
    }
    
}

::SUNMatrix SUNMatrixWrapper::SUNMatrix() const {
    return matrix;
}
} // namespace amici
