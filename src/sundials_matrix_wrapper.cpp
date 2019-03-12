#include <amici/sundials_matrix_wrapper.h>

#include <amici/cblas.h>

#include <new> // bad_alloc
#include <utility>
#include <stdexcept> // invalid_argument and domain_error

namespace amici {

SUNMatrixWrapper::SUNMatrixWrapper(int M, int N, int NNZ, int sparsetype)
    : matrix(SUNSparseMatrix(M, N, NNZ, sparsetype)) {
    if (sparsetype != CSC_MAT && sparsetype != CSR_MAT)
        throw std::invalid_argument("Invalid sparsetype. Must be CSC_MAT or "
                                    "CSR_MAT");

    if (NNZ && !matrix)
        throw std::bad_alloc();
    
    if (NNZ) {
        data_ptr = SM_DATA_S(matrix);
        indexptrs_ptr = SM_INDEXPTRS_S(matrix);
        indexvals_ptr = SM_INDEXVALS_S(matrix);
    }
}

SUNMatrixWrapper::SUNMatrixWrapper(int M, int N)
    : matrix(SUNDenseMatrix(M, N)) {
    if (M && N && !matrix)
        throw std::bad_alloc();
    if (M && N)
        data_ptr = SM_DATA_D(matrix);
}

SUNMatrixWrapper::SUNMatrixWrapper(int M, int ubw, int lbw)
    : matrix(SUNBandMatrix(M, ubw, lbw)) {
    if (M && !matrix)
        throw std::bad_alloc();
    if (M)
        data_ptr = SM_DATA_B(matrix);
}

SUNMatrixWrapper::SUNMatrixWrapper(const SUNMatrixWrapper &A, realtype droptol,
                                   int sparsetype) {
    if (sparsetype != CSC_MAT && sparsetype != CSR_MAT)
        throw std::invalid_argument("Invalid sparsetype. Must be CSC_MAT or "
                                    "CSR_MAT");

    switch (SUNMatGetID(A.get())) {
    case SUNMATRIX_DENSE:
        matrix = SUNSparseFromDenseMatrix(A.get(), droptol, sparsetype);
        break;
    case SUNMATRIX_BAND:
        matrix = SUNSparseFromBandMatrix(A.get(), droptol, sparsetype);
        break;
    default:
        throw std::invalid_argument("Invalid Matrix. Must be SUNMATRIX_DENSE or"
                                    " SUNMATRIX_BAND");
    }

    if (!matrix)
        throw std::bad_alloc();
    else {
        data_ptr = SM_DATA_S(matrix);
        indexptrs_ptr = SM_INDEXPTRS_S(matrix);
        indexvals_ptr = SM_INDEXVALS_S(matrix);
    }
}

SUNMatrixWrapper::SUNMatrixWrapper(SUNMatrix mat) : matrix(mat) {}

SUNMatrixWrapper::~SUNMatrixWrapper() {
    if (matrix)
        SUNMatDestroy(matrix);
}

SUNMatrixWrapper::SUNMatrixWrapper(const SUNMatrixWrapper &other) {
    if (!other.matrix)
        return;

    matrix = SUNMatClone(other.matrix);
    if (!matrix)
        throw std::bad_alloc();

    SUNMatCopy(other.matrix, matrix);
    data_ptr = other.data_ptr;
    indexptrs_ptr = other.indexptrs_ptr;
    indexvals_ptr = other.indexvals_ptr;
}

SUNMatrixWrapper::SUNMatrixWrapper(SUNMatrixWrapper &&other) noexcept {
    std::swap(matrix, other.matrix);
    data_ptr = other.data_ptr;
    indexptrs_ptr = other.indexptrs_ptr;
    indexvals_ptr = other.indexvals_ptr;
}

SUNMatrixWrapper &SUNMatrixWrapper::operator=(const SUNMatrixWrapper &other) {
    return *this = SUNMatrixWrapper(other);
}

SUNMatrixWrapper &SUNMatrixWrapper::
operator=(SUNMatrixWrapper &&other) noexcept {
    std::swap(matrix, other.matrix);
    data_ptr = other.data_ptr;
    indexptrs_ptr = other.indexptrs_ptr;
    indexvals_ptr = other.indexvals_ptr;
    return *this;
}

realtype *SUNMatrixWrapper::data() const {
    return data_ptr;
}

sunindextype SUNMatrixWrapper::rows() const {
    if (!matrix)
        return 0;

    switch (SUNMatGetID(matrix)) {
    case SUNMATRIX_DENSE:
        return SM_ROWS_D(matrix);
    case SUNMATRIX_SPARSE:
        return SM_ROWS_S(matrix);
    case SUNMATRIX_BAND:
        return SM_ROWS_B(matrix);
    case SUNMATRIX_CUSTOM:
        throw std::domain_error("Amici currently does not support custom matrix"
                                " types.");
    default:
        throw std::domain_error("Invalid SUNMatrix type.");
    }
}

sunindextype SUNMatrixWrapper::columns() const {
    if (!matrix)
        return 0;

    switch (SUNMatGetID(matrix)) {
    case SUNMATRIX_DENSE:
        return SM_COLUMNS_D(matrix);
    case SUNMATRIX_SPARSE:
        return SM_COLUMNS_S(matrix);
    case SUNMATRIX_BAND:
        return SM_COLUMNS_B(matrix);
    case SUNMATRIX_CUSTOM:
        throw std::domain_error("Amici currently does not support custom matrix"
                                " types.");
    default:
        throw std::domain_error("Invalid SUNMatrix type.");
    }
}

sunindextype *SUNMatrixWrapper::indexvals() const {
    return indexvals_ptr;
}

sunindextype *SUNMatrixWrapper::indexptrs() const {
    return indexptrs_ptr;
}

int SUNMatrixWrapper::sparsetype() const {
    if (SUNMatGetID(matrix) == SUNMATRIX_SPARSE)
        return SM_SPARSETYPE_S(matrix);
    else
        throw std::domain_error("Function only available for sparse matrices");
}

void SUNMatrixWrapper::reset() {
    if (matrix)
        SUNMatZero(matrix);
}

void SUNMatrixWrapper::multiply(std::vector<realtype> &c,
                                const std::vector<realtype> &b) const {
    if (static_cast<sunindextype>(c.size()) != rows())
        throw std::invalid_argument("Dimension mismatch between number of rows"
                                    "in A and elements in c");

    if (static_cast<sunindextype>(b.size()) != columns())
        throw std::invalid_argument("Dimension mismatch between number of cols"
                                    "in A and elements in b");

    multiply(c.data(), b.data());
}

void SUNMatrixWrapper::multiply(N_Vector c, const N_Vector b) const {
    if (NV_LENGTH_S(c) != rows())
        throw std::invalid_argument("Dimension mismatch between number of rows"
                                    "in A and elements in c");

    if (NV_LENGTH_S(b) != columns())
        throw std::invalid_argument("Dimension mismatch between number of cols"
                                    "in A and elements in b");

    multiply(NV_DATA_S(c), NV_DATA_S(b));
}

void SUNMatrixWrapper::multiply(realtype *c, const realtype *b) const {
    if (!matrix)
        return;

    switch (SUNMatGetID(matrix)) {
    case SUNMATRIX_DENSE:
        amici_dgemv(BLASLayout::colMajor, BLASTranspose::noTrans, rows(),
                    columns(), 1.0, data(), rows(), b, 1, 1.0, c, 1);
        break;
    case SUNMATRIX_SPARSE:

        switch (sparsetype()) {
        case CSC_MAT:
            for (sunindextype i = 0; i < columns(); ++i) {
                for (sunindextype k = indexptrs()[i]; k < indexptrs()[i + 1];
                     ++k) {
                    c[indexvals()[k]] += data()[k] * b[i];
                }
            }
            break;
        case CSR_MAT:
            for (sunindextype i = 0; i < rows(); ++i) {
                for (sunindextype k = indexptrs()[i]; k < indexptrs()[i + 1];
                     ++k) {
                    c[i] += data()[k] * b[indexvals()[k]];
                }
            }
            break;
        }
        break;
    case SUNMATRIX_BAND:
        throw std::domain_error("Not Implemented.");
    case SUNMATRIX_CUSTOM:
        throw std::domain_error("Amici currently does not support custom"
                                " matrix types.");
    }
}

SUNMatrix SUNMatrixWrapper::get() const { return matrix; }

} // namespace amici

