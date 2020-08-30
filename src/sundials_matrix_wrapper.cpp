#include <amici/sundials_matrix_wrapper.h>

#include <amici/cblas.h>

#include <new> // bad_alloc
#include <utility>
#include <stdexcept> // invalid_argument and domain_error

namespace amici {

SUNMatrixWrapper::SUNMatrixWrapper(int M, int N, int NNZ, int sparsetype)
    : matrix_(SUNSparseMatrix(M, N, NNZ, sparsetype)) {

    if (sparsetype != CSC_MAT && sparsetype != CSR_MAT)
        throw std::invalid_argument("Invalid sparsetype. Must be CSC_MAT or "
                                    "CSR_MAT");

    if (NNZ && !matrix_)
        throw std::bad_alloc();

    update_ptrs();
}

SUNMatrixWrapper::SUNMatrixWrapper(int M, int N)
    : matrix_(SUNDenseMatrix(M, N)) {
    if (M && N && !matrix_)
        throw std::bad_alloc();

    update_ptrs();
}

SUNMatrixWrapper::SUNMatrixWrapper(int M, int ubw, int lbw)
    : matrix_(SUNBandMatrix(M, ubw, lbw)) {
    if (M && !matrix_)
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
        matrix_ = SUNSparseFromDenseMatrix(A.get(), droptol, sparsetype);
        break;
    case SUNMATRIX_BAND:
        matrix_ = SUNSparseFromBandMatrix(A.get(), droptol, sparsetype);
        break;
    default:
        throw std::invalid_argument("Invalid Matrix. Must be SUNMATRIX_DENSE or"
                                    " SUNMATRIX_BAND");
    }

    if (!matrix_)
        throw std::bad_alloc();

    update_ptrs();
}

SUNMatrixWrapper::SUNMatrixWrapper(SUNMatrix mat) : matrix_(mat) {
    update_ptrs();
}

SUNMatrixWrapper::~SUNMatrixWrapper() {
    if (matrix_)
        SUNMatDestroy(matrix_);
}

SUNMatrixWrapper::SUNMatrixWrapper(const SUNMatrixWrapper &other) {
    if (!other.matrix_)
        return;

    matrix_ = SUNMatClone(other.matrix_);
    if (!matrix_)
        throw std::bad_alloc();

    SUNMatCopy(other.matrix_, matrix_);
    update_ptrs();
}

SUNMatrixWrapper::SUNMatrixWrapper(SUNMatrixWrapper &&other) {
    std::swap(matrix_, other.matrix_);
    update_ptrs();
}

SUNMatrixWrapper &SUNMatrixWrapper::operator=(const SUNMatrixWrapper &other) {
    if(&other == this)
        return *this;
    return *this = SUNMatrixWrapper(other);
}

SUNMatrixWrapper &SUNMatrixWrapper::
operator=(SUNMatrixWrapper &&other) {
    std::swap(matrix_, other.matrix_);
    update_ptrs();
    return *this;
}

void SUNMatrixWrapper::reallocate(int NNZ) {
    if (sparsetype() != CSC_MAT && sparsetype() != CSR_MAT)
        throw std::invalid_argument("Invalid sparsetype. Must be CSC_MAT or "
                                    "CSR_MAT");
    
    if (SUNSparseMatrix_Reallocate(matrix_, NNZ) != AMICI_SUCCESS)
        throw std::bad_alloc();
    
    if ((NNZ && !matrix_) | (NNZ != nonzeros()))
        throw std::bad_alloc();
    update_ptrs();
}

void SUNMatrixWrapper::realloc() {
    if (sparsetype() != CSC_MAT && sparsetype() != CSR_MAT)
        throw std::invalid_argument("Invalid sparsetype. Must be CSC_MAT or "
                                    "CSR_MAT");
    if (SUNSparseMatrix_Realloc(matrix_) != AMICI_SUCCESS)
        throw std::bad_alloc();
        
    if (nonzeros() && !matrix_)
        throw std::bad_alloc();
    update_ptrs();
}

realtype *SUNMatrixWrapper::data() const {
    return data_ptr_;
}

sunindextype SUNMatrixWrapper::rows() const {
    if (!matrix_)
        return 0;

    switch (SUNMatGetID(matrix_)) {
    case SUNMATRIX_DENSE:
        return SM_ROWS_D(matrix_);
    case SUNMATRIX_SPARSE:
        return SM_ROWS_S(matrix_);
    case SUNMATRIX_BAND:
        return SM_ROWS_B(matrix_);
    case SUNMATRIX_CUSTOM:
        throw std::domain_error("Amici currently does not support custom matrix"
                                " types.");
    default:
        throw std::domain_error("Invalid SUNMatrix type.");
    }
}

sunindextype SUNMatrixWrapper::columns() const {
    if (!matrix_)
        return 0;

    switch (SUNMatGetID(matrix_)) {
    case SUNMATRIX_DENSE:
        return SM_COLUMNS_D(matrix_);
    case SUNMATRIX_SPARSE:
        return SM_COLUMNS_S(matrix_);
    case SUNMATRIX_BAND:
        return SM_COLUMNS_B(matrix_);
    case SUNMATRIX_CUSTOM:
        throw std::domain_error("Amici currently does not support custom matrix"
                                " types.");
    default:
        throw std::domain_error("Invalid SUNMatrix type.");
    }
}

sunindextype SUNMatrixWrapper::nonzeros() const {
    if (!matrix_)
        return 0;

    switch (SUNMatGetID(matrix_)) {
    case SUNMATRIX_SPARSE:
        return SM_NNZ_S(matrix_);
    default:
        throw std::domain_error("Non-zeros property only available for "
                                "sparse matrices");
    }
}

sunindextype *SUNMatrixWrapper::indexvals() const {
    return indexvals_ptr_;
}

sunindextype *SUNMatrixWrapper::indexptrs() const {
    return indexptrs_ptr_;
}

int SUNMatrixWrapper::sparsetype() const {
    if (SUNMatGetID(matrix_) == SUNMATRIX_SPARSE)
        return SM_SPARSETYPE_S(matrix_);
    throw std::domain_error("Function only available for sparse matrices");
}

void SUNMatrixWrapper::reset() {
    if (matrix_)
        SUNMatZero(matrix_);
}

void SUNMatrixWrapper::scale(realtype a) {
    if (matrix_) {
        int nonzeros_ = static_cast<int>(nonzeros());
        for (int i = 0; i < nonzeros_; ++i)
            data_ptr_[i] *= a;
    }
}

void SUNMatrixWrapper::multiply(N_Vector c, const_N_Vector b) const {
    multiply(gsl::make_span<realtype>(NV_DATA_S(c), NV_LENGTH_S(c)),
             gsl::make_span<const realtype>(NV_DATA_S(b), NV_LENGTH_S(b)));
}

void SUNMatrixWrapper::multiply(gsl::span<realtype> c,
                                gsl::span<const realtype> b) const {
    if (!matrix_)
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

    switch (SUNMatGetID(matrix_)) {
    case SUNMATRIX_DENSE:
        amici_dgemv(BLASLayout::colMajor, BLASTranspose::noTrans,
                    static_cast<int>(nrows), static_cast<int>(ncols),
                    1.0, data(), static_cast<int>(nrows),
                    b.data(), 1, 1.0, c.data(), 1);
        break;
    case SUNMATRIX_SPARSE:
        if(!SM_NNZ_S(matrix_)) {
            /* empty matrix, nothing to multiply, return to avoid out-of-bounds
             * access of pointer access below
             */
            return;
        }
        switch (sparsetype()) {
        case CSC_MAT:
            for (sunindextype i = 0; i < ncols; ++i) {
                for (sunindextype k = indexptrs_ptr_[i]; k < indexptrs_ptr_[i + 1];
                     ++k) {
                    c[indexvals_ptr_[k]] += data_ptr_[k] * b[i];
                }
            }
            break;
        case CSR_MAT:
            for (sunindextype i = 0; i < nrows; ++i) {
                for (sunindextype k = indexptrs_ptr_[i]; k < indexptrs_ptr_[i + 1];
                     ++k) {
                    c[i] += data_ptr_[k] * b[indexvals_ptr_[k]];
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
    case SUNMATRIX_SLUNRLOC:
        throw std::domain_error("Not Implemented.");
    case SUNMATRIX_CUSPARSE:
        throw std::domain_error("Not Implemented.");
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
    if (!matrix_)
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

    if (SUNMatGetID(matrix_) != SUNMATRIX_SPARSE)
        throw std::invalid_argument("Reordered multiply only implemented for "
                                    "sparse matrices, but A is not sparse");

    if (sparsetype() != CSC_MAT)
        throw std::invalid_argument("Reordered multiply only implemented for "
                                    "matrix type CSC, but A is not of type CSC");

    /* Carry out actual multiplication */
    if (transpose) {
        for (int i = 0; i < (int)cols.size(); ++i)
            for (sunindextype k = indexptrs_ptr_[cols[i]];
                 k < indexptrs_ptr_[cols[i] + 1]; ++k)
                c[i] += data_ptr_[k] * b[indexvals_ptr_[k]];
    } else {
        for (sunindextype i = 0; i < ncols; ++i)
            for (sunindextype k = indexptrs_ptr_[cols[i]];
                 k < indexptrs_ptr_[cols[i] + 1]; ++k)
                c[indexvals_ptr_[k]] += data_ptr_[k] * b[i];
    }
}


void SUNMatrixWrapper::sparse_multiply(SUNMatrixWrapper *C,
                                       SUNMatrixWrapper *B) const {
    if (!matrix_)
        return;

    sunindextype nrows = rows();
    sunindextype ncols = columns();

    if (SUNMatGetID(matrix_) != SUNMATRIX_SPARSE)
        throw std::invalid_argument("Matrix A not sparse in sparse_multiply");

    if (sparsetype() != CSC_MAT)
        throw std::invalid_argument("Matrix A not of type CSC_MAT");

    if (SUNMatGetID(B->matrix_) != SUNMATRIX_SPARSE)
        throw std::invalid_argument("Matrix B not sparse in sparse_multiply");

    if (B->sparsetype() != CSC_MAT)
        throw std::invalid_argument("Matrix B not of type CSC_MAT");

    if (SUNMatGetID(C->matrix_) != SUNMATRIX_SPARSE)
        throw std::invalid_argument("Matrix C not sparse in sparse_multiply");

    if (C->sparsetype() != CSC_MAT)
        throw std::invalid_argument("Matrix C not of type CSC_MAT");

    if (C->rows() != nrows)
        throw std::invalid_argument("Dimension mismatch between number of rows "
                                    "in A (" + std::to_string(nrows)
                                    + ") and number of rows in C ("
                                    + std::to_string((int)C->rows()) + ")");
    
    if (C->columns() != B->columns())
        throw std::invalid_argument("Dimension mismatch between number of cols "
                                    "in B (" + std::to_string(B->columns())
                                    + ") and  number of cols in C ("
                                    + std::to_string((int)C->columns()) + ")");

    if (B->rows() != ncols)
        throw std::invalid_argument("Dimension mismatch between number of cols "
                                    "in A (" + std::to_string(ncols)
                                    + ") and number of rows in B ("
                                    + std::to_string((int)B->rows()) + ")");
    
    if (ncols == 0 || nrows == 0 || B->columns() == 0)
        return; // matrix will also have zero size
    
    if (nonzeros() == 0 || B->nonzeros() == 0)
        return; //  nothing to multiply
    

    /* see https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/master/CSparse/Source/cs_multiply.c
     * modified such that we don't need to use CSparse memory structure and can
     * work with preallocated C. This should minimize number of necessary
     * reallocations as we can assume that C doesn't change size.
     */
    
    sunindextype nz = 0; // this keeps track of the nonzero index in C
    auto Bx = B->data();
    auto Bi = B->indexvals();
    auto Bp = B->indexptrs();
    
    sunindextype j, p;
    auto Cx = C->data();
    auto Ci = C->indexvals();
    auto Cp = C->indexptrs();
    
    sunindextype m = nrows;
    sunindextype n = B->columns();
    
    auto w = std::vector<sunindextype>(m);
    auto x = std::vector<realtype>(m);
    
    for (j = 0; j < n; j++)
    {
        Cp[j] = nz;                          /* column j of C starts here */
        // CSparse reallocation happend when it _might_ have been necessary,
        // which could happe in every call to mulitply_sparse, even though
        // allocated space was actually sufficient.
        for (p = Bp[j]; p < Bp[j+1]; p++) {
            nz = scatter(Bi[p], Bx[p], w.data(), x.data(), j+1, C, nz);
            // in case of reallocation of C in scatter
            Cx = C->data(); Ci = C->indexvals(); Cp = C->indexptrs();
        }
        for (p = Cp[j]; p < nz; p++)
            Cx[p] = x.at(Ci[p]); // copy data to C
    }
    Cp[n] = nz;
    C->realloc();
}

sunindextype SUNMatrixWrapper::scatter(const sunindextype j,
                                       const realtype beta,
                                       sunindextype *w, realtype *x,
                                       const sunindextype mark,
                                       SUNMatrixWrapper *C,
                                       sunindextype nz) const {
    if (sparsetype() != CSC_MAT)
        throw std::invalid_argument("Matrix A not of type CSC_MAT");
    
    if (C->sparsetype() != CSC_MAT)
        throw std::invalid_argument("Matrix C not of type CSC_MAT");
    
    /* see https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/master/CSparse/Source/cs_scatter.c */
    
    auto Ci = C->indexvals();
    auto Ap = indexptrs();
    auto Ai = indexvals();
    auto Ax = data();
    for (sunindextype p = Ap[j]; p < Ap[j+1]; p++)
    {
        auto i = Ai[p];                   /* A(i,j) is nonzero */
        if (w[i] < mark) {
            if (nz + 1 > C->nonzeros())
            {
                // move reallocation here since we can explicitely only call
                // this when space in C is actually not sufficient.
                C->reallocate(2*C->nonzeros()+1);
                Ci = C->indexvals();
            }
            w[i] = mark;                  /* i is new entry in column j */
            Ci[nz++] = i;                 /* add i to pattern of C(:,j) */
            x[i] = beta * Ax[p];          /* x(i) = beta*A(i,j) */
        }
        else x[i] += beta * Ax[p];        /* i exists in C(:,j) already */
    }
    return nz;
}

void SUNMatrixWrapper::zero()
{
    if(int res = SUNMatZero(matrix_))
        throw std::runtime_error("SUNMatrixWrapper::zero() failed with "
                                 + std::to_string(res));
}

void SUNMatrixWrapper::update_ptrs() {
    if(!matrix_)
        return;

    switch (SUNMatGetID(matrix_)) {
    case SUNMATRIX_DENSE:
        if (columns() > 0 && rows() > 0)
            data_ptr_ = SM_DATA_D(matrix_);
        break;
    case SUNMATRIX_SPARSE:
        if (SM_NNZ_S(matrix_) > 0) {
            data_ptr_ = SM_DATA_S(matrix_);
            indexptrs_ptr_ = SM_INDEXPTRS_S(matrix_);
            indexvals_ptr_ = SM_INDEXVALS_S(matrix_);
        }
        break;
    case SUNMATRIX_BAND:
        if (columns() > 0 && rows() > 0)
            data_ptr_ = SM_DATA_B(matrix_);
        break;
    case SUNMATRIX_CUSTOM:
        throw std::domain_error("Amici currently does not support "
                                "custom matrix types.");
    case SUNMATRIX_SLUNRLOC:
        throw std::domain_error("Not Implemented.");
    case SUNMATRIX_CUSPARSE:
        throw std::domain_error("Not Implemented.");
    }
}

SUNMatrix SUNMatrixWrapper::get() const { return matrix_; }

} // namespace amici

