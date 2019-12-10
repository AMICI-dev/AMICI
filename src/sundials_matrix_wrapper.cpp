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

    update_ptrs();
}

SUNMatrixWrapper::SUNMatrixWrapper(int M, int N)
    : matrix(SUNDenseMatrix(M, N)) {
    if (M && N && !matrix)
        throw std::bad_alloc();

    update_ptrs();
}

SUNMatrixWrapper::SUNMatrixWrapper(int M, int ubw, int lbw)
    : matrix(SUNBandMatrix(M, ubw, lbw)) {
    if (M && !matrix)
        throw std::bad_alloc();

    update_ptrs();
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

    update_ptrs();
}

SUNMatrixWrapper::SUNMatrixWrapper(SUNMatrix mat) : matrix(mat) {
    update_ptrs();
}

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
    update_ptrs();
}

SUNMatrixWrapper::SUNMatrixWrapper(SUNMatrixWrapper &&other) {
    std::swap(matrix, other.matrix);
    update_ptrs();
}

SUNMatrixWrapper &SUNMatrixWrapper::operator=(const SUNMatrixWrapper &other) {
    if(&other == this)
        return *this;
    return *this = SUNMatrixWrapper(other);
}

SUNMatrixWrapper &SUNMatrixWrapper::
operator=(SUNMatrixWrapper &&other) {
    std::swap(matrix, other.matrix);
    update_ptrs();
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

sunindextype SUNMatrixWrapper::nonzeros() const {
    if (!matrix)
        return 0;

    switch (SUNMatGetID(matrix)) {
    case SUNMATRIX_SPARSE:
        return SM_NNZ_S(matrix);
    default:
        throw std::domain_error("Non-zeros property only available for "
                                "sparse matrices");
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
    throw std::domain_error("Function only available for sparse matrices");
}

void SUNMatrixWrapper::reset() {
    if (matrix)
        SUNMatZero(matrix);
}

void SUNMatrixWrapper::multiply(N_Vector c, const_N_Vector b) const {
    multiply(gsl::make_span<realtype>(NV_DATA_S(c), NV_LENGTH_S(c)),
             gsl::make_span<const realtype>(NV_DATA_S(b), NV_LENGTH_S(b)));
}

void SUNMatrixWrapper::multiply(gsl::span<realtype> c,
                                gsl::span<const realtype> b) const {
    if (!matrix)
        return;

    sunindextype nrows = rows();
    sunindextype ncols = columns();

    if (static_cast<sunindextype>(c.size()) != nrows)
        throw std::invalid_argument("Dimension mismatch between number of rows "
                                    "in A (" + std::to_string(nrows) + ") and "
                                    "elements in c (" + std::to_string(c.size())
                                    + ")");

    if (static_cast<sunindextype>(b.size()) != ncols)
        throw std::invalid_argument("Dimension mismatch between number of cols "
                                    "in A (" + std::to_string(ncols)
                                    + ") and elements in b ("
                                    + std::to_string(b.size()) + ")");

    switch (SUNMatGetID(matrix)) {
    case SUNMATRIX_DENSE:
        amici_dgemv(BLASLayout::colMajor, BLASTranspose::noTrans, nrows,
                    ncols, 1.0, data(), nrows, b.data(), 1, 1.0, c.data(), 1);
        break;
    case SUNMATRIX_SPARSE:

        switch (sparsetype()) {
        case CSC_MAT:
            for (sunindextype i = 0; i < ncols; ++i) {
                for (sunindextype k = indexptrs_ptr[i]; k < indexptrs_ptr[i + 1];
                     ++k) {
                    c[indexvals_ptr[k]] += data_ptr[k] * b[i];
                }
            }
            break;
        case CSR_MAT:
            for (sunindextype i = 0; i < nrows; ++i) {
                for (sunindextype k = indexptrs_ptr[i]; k < indexptrs_ptr[i + 1];
                     ++k) {
                    c[i] += data_ptr[k] * b[indexvals_ptr[k]];
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

void SUNMatrixWrapper::multiply(N_Vector c,
                                const N_Vector b,
                                gsl::span <const int> cols,
                                bool transpose) const {
    multiply(gsl::make_span<realtype>(NV_DATA_S(c), NV_LENGTH_S(c)),
             gsl::make_span<const realtype>(NV_DATA_S(b), NV_LENGTH_S(b)),
             cols, transpose);
}

void SUNMatrixWrapper::multiply(gsl::span<realtype> c,
                                gsl::span<const realtype> b,
                                gsl::span<const int> cols,
                                bool transpose) const {
    if (!matrix)
        return;

    sunindextype nrows = rows();
    sunindextype ncols = columns();

    if (transpose) {
        if (c.size() != cols.size())
            throw std::invalid_argument("Dimension mismatch between number of cols "
                                        "in index vector cols (" + std::to_string(ncols)
                                        + ") and elements in c (" + std::to_string(c.size())
                                        + " ), when using transposed A");

        if (static_cast<sunindextype>(b.size()) != nrows)
            throw std::invalid_argument("Dimension mismatch between number of rows "
                                        "in A (" + std::to_string(nrows) + ") and "
                                        "elements in b (" + std::to_string(b.size())
                                        + "), when using transposed A");
    } else {
        if (static_cast<sunindextype>(c.size()) != nrows)
            throw std::invalid_argument("Dimension mismatch between number of rows "
                                        "in A (" + std::to_string(nrows) + ") and "
                                        "elements in c (" + std::to_string(c.size())
                                        + ")");

        if (static_cast<sunindextype>(b.size()) != ncols)
            throw std::invalid_argument("Dimension mismatch between number of cols "
                                        "in A (" + std::to_string(ncols)
                                        + ") and elements in b ("
                                        + std::to_string(b.size()) + ")");
    }

    if (SUNMatGetID(matrix) != SUNMATRIX_SPARSE)
        throw std::invalid_argument("Reordered multiply only implemented for "
                                    "sparse matrices, but A is not sparse");

    if (sparsetype() != CSC_MAT)
        throw std::invalid_argument("Reordered multiply only implemented for "
                                    "matrix type CSC, but A is not of type CSC");

    /* Carry out actual multiplication */
    if (transpose) {
        for (int i = 0; i < (int)cols.size(); ++i)
            for (sunindextype k = indexptrs_ptr[cols[i]];
                 k < indexptrs_ptr[cols[i] + 1]; ++k)
                c[i] += data_ptr[k] * b[indexvals_ptr[k]];
    } else {
        for (sunindextype i = 0; i < ncols; ++i)
            for (sunindextype k = indexptrs_ptr[cols[i]];
                 k < indexptrs_ptr[cols[i] + 1]; ++k)
                c[indexvals_ptr[k]] += data_ptr[k] * b[i];
    }
}


void SUNMatrixWrapper::sparse_multiply(SUNMatrixWrapper *C,
                                       SUNMatrixWrapper *B) const {
    if (!matrix)
        return;

    sunindextype nrows = rows();
    sunindextype ncols = columns();

    if (SUNMatGetID(matrix) != SUNMATRIX_SPARSE)
        throw std::invalid_argument("Matrix A not sparse in sparse_multiply");

    if (sparsetype() != CSC_MAT)
        throw std::invalid_argument("Matrix A not of type CSC_MAT");

    if (SUNMatGetID(B->matrix) != SUNMATRIX_SPARSE)
        throw std::invalid_argument("Matrix B not sparse in sparse_multiply");

    if (B->sparsetype() != CSC_MAT)
        throw std::invalid_argument("Matrix B not of type CSC_MAT");

    if (SUNMatGetID(C->matrix) != SUNMATRIX_SPARSE)
        throw std::invalid_argument("Matrix C not sparse in sparse_multiply");

    if (C->sparsetype() != CSC_MAT)
        throw std::invalid_argument("Matrix C not of type CSC_MAT");

    if (C->rows() != nrows)
        throw std::invalid_argument("Dimension mismatch between number of rows "
                                    "in A (" + std::to_string(nrows) + ") and "
                                    "number of rows in C ("
                                    + std::to_string((int)C->rows()) + ")");

    if (B->rows() != ncols)
        throw std::invalid_argument("Dimension mismatch between number of rows "
                                    "in A (" + std::to_string(ncols)
                                    + ") and number of cols in B ("
                                    + std::to_string((int)B->rows()) + ")");

    /* Carry out actual multiplication */
    sunindextype iC_data = 0;
    for (sunindextype iC_col = 0; iC_col < C->columns(); ++iC_col) {
        for(sunindextype iC_row = C->indexptrs_ptr[iC_col];
            iC_row < C->indexptrs_ptr[iC_col + 1]; ++iC_row) {
            // Current entry in C: (C.indexvals[iC_row], iC_col)
            for(sunindextype iB_row = B->indexptrs_ptr[iC_col];
                iB_row < B->indexptrs_ptr[iC_col + 1]; ++iB_row) {
                // Loop over column iC_col in B
                sunindextype iA_col = B->indexvals_ptr[iB_row];
                for (sunindextype iA_row = indexptrs_ptr[iA_col];
                     iA_row < indexptrs_ptr[iA_col + 1]; ++iA_row)
                    // loop over entries in column col_in_A
                    // if two entries match together, carry out multiplication
                    if (indexvals_ptr[iA_row] == C->indexvals_ptr[iC_row])
                        C->data_ptr[iC_data] +=
                            data_ptr[iA_row] * B->data_ptr[iB_row];
            }
            iC_data++;
        }
    }
}

void SUNMatrixWrapper::zero()
{
    if(int res = SUNMatZero(matrix))
        throw std::runtime_error("SUNMatrixWrapper::zero() failed with "
                                 + std::to_string(res));
}

void SUNMatrixWrapper::update_ptrs() {
    if(!matrix)
        return;

    switch (SUNMatGetID(matrix)) {
    case SUNMATRIX_DENSE:
        if (columns() > 0 && rows() > 0)
            data_ptr = SM_DATA_D(matrix);
        break;
    case SUNMATRIX_SPARSE:
        if (SM_NNZ_S(matrix) > 0) {
            data_ptr = SM_DATA_S(matrix);
            indexptrs_ptr = SM_INDEXPTRS_S(matrix);
            indexvals_ptr = SM_INDEXVALS_S(matrix);
        }
        break;
    case SUNMATRIX_BAND:
        if (columns() > 0 && rows() > 0)
            data_ptr = SM_DATA_B(matrix);
        break;
    case SUNMATRIX_CUSTOM:
        throw std::domain_error("Amici currently does not support "
                                "custom matrix types.");
    }
}

SUNMatrix SUNMatrixWrapper::get() const { return matrix; }

} // namespace amici

