#include <amici/sundials_matrix_wrapper.h>

#include <utility>
#include <new> // bad_alloc

#include <sundials/sundials_dense.h> // DenseCopy

namespace amici {

SlsMatWrapper::SlsMatWrapper(int M, int N, int NNZ, int sparsetype)
    : matrix(SparseNewMat(M, N, NNZ, sparsetype))
{
    if(NNZ && !matrix)
        throw std::bad_alloc();
}

SlsMatWrapper::SlsMatWrapper(SlsMat mat)
    : matrix(mat)
{}

SlsMatWrapper::~SlsMatWrapper()
{
    if(matrix)
        SparseDestroyMat(matrix);
}

SlsMatWrapper::SlsMatWrapper(const SlsMatWrapper &other)
{
    if(!other.matrix)
        return;

    matrix = SparseNewMat(other.matrix->M,
                          other.matrix->N,
                          other.matrix->NNZ,
                          other.matrix->sparsetype);
    if(!matrix)
        throw std::bad_alloc();

    SparseCopyMat(other.matrix, matrix);
}

SlsMatWrapper::SlsMatWrapper(SlsMatWrapper &&other) noexcept
{
    std::swap(matrix, other.matrix);

}

SlsMatWrapper &SlsMatWrapper::operator=(const SlsMatWrapper &other)
{
    return *this = SlsMatWrapper(other);
}

SlsMatWrapper &SlsMatWrapper::operator=(SlsMatWrapper &&other) noexcept
{
    std::swap(matrix, other.matrix);
    return *this;
}

realtype *SlsMatWrapper::data() {
    if(matrix)
        return matrix->data;
    return nullptr;
}

SlsMat SlsMatWrapper::slsmat() const {
    return matrix;
}



DlsMatWrapper::DlsMatWrapper(long int M, long int N)
    : matrix(NewDenseMat(M, N))
{
    if((M*N > 0) && !matrix)
        throw std::bad_alloc();
}

DlsMatWrapper::DlsMatWrapper(DlsMat mat)
    : matrix(mat)
{}

DlsMatWrapper::~DlsMatWrapper()
{
    if(matrix)
        DestroyMat(matrix);
}

DlsMatWrapper::DlsMatWrapper(const DlsMatWrapper &other)
{
    if(!other.matrix)
        return;

    matrix = NewDenseMat(other.matrix->M,
                          other.matrix->N);
    if(!matrix)
        throw std::bad_alloc();

    DenseCopy(other.matrix, matrix);
}

DlsMatWrapper::DlsMatWrapper(DlsMatWrapper &&other) noexcept
{
    std::swap(matrix, other.matrix);
}

DlsMatWrapper &DlsMatWrapper::operator=(const DlsMatWrapper &other)
{
    return *this = DlsMatWrapper(other);
}

DlsMatWrapper &DlsMatWrapper::operator=(DlsMatWrapper &&other) noexcept
{
    std::swap(matrix, other.matrix);
    return *this;
}

realtype *DlsMatWrapper::data() {
    if(matrix)
        return matrix->data;
    return nullptr;
}

DlsMat DlsMatWrapper::dlsmat() const {
    return matrix;
}
} // namespace amici
