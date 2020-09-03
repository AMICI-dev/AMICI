#include <amici/sundials_matrix_wrapper.h>
#include <sundials/sundials_matrix.h> // return codes

#include <amici/cblas.h>

#include <new> // bad_alloc
#include <utility>
#include <stdexcept> // invalid_argument and domain_error
#include <assert.h>

namespace amici {

SUNMatrixWrapper::SUNMatrixWrapper(int M, int N, int NNZ, int sparsetype)
    : matrix_(SUNSparseMatrix(M, N, NNZ, sparsetype)) {

    if (sparsetype != CSC_MAT && sparsetype != CSR_MAT)
        throw std::invalid_argument("Invalid sparsetype. Must be CSC_MAT or "
                                    "CSR_MAT");

    if (NNZ && M && N && !matrix_)
        throw std::bad_alloc();
        
    assert(num_nonzeros() == 0);
    assert(NNZ == capacity());
    assert(M == rows() || !matrix_);
    assert(N == columns() || !matrix_);

    update_ptrs();
}

SUNMatrixWrapper::SUNMatrixWrapper(int M, int N)
    : matrix_(SUNDenseMatrix(M, N)) {
    if (M && N && !matrix_)
        throw std::bad_alloc();
        
    assert(M == rows());
    assert(N == columns());

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

void SUNMatrixWrapper::reallocate(sunindextype NNZ) {
    if (sparsetype() != CSC_MAT && sparsetype() != CSR_MAT)
        throw std::invalid_argument("Invalid sparsetype. Must be CSC_MAT or "
                                    "CSR_MAT.");
    
    if (int ret = SUNSparseMatrix_Reallocate(matrix_, NNZ) != SUNMAT_SUCCESS)
        throw std::runtime_error("SUNSparseMatrix_Reallocate failed with "
                                 "error code " + std::to_string(ret) + ".");

    update_ptrs();
    assert((NNZ && columns() && rows()) ^ !matrix_);
    assert(NNZ == capacity());
}

void SUNMatrixWrapper::realloc() {
    if (sparsetype() != CSC_MAT && sparsetype() != CSR_MAT)
        throw std::invalid_argument("Invalid sparsetype. Must be CSC_MAT or "
                                    "CSR_MAT.");
    if (int ret = SUNSparseMatrix_Realloc(matrix_) != SUNMAT_SUCCESS)
        throw std::runtime_error("SUNSparseMatrix_Realloc failed with "
                                 "error code " + std::to_string(ret) + ".");
    
    update_ptrs();
    assert(capacity() ^ !matrix_);
    assert(capacity() == num_nonzeros());
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

sunindextype SUNMatrixWrapper::capacity() const {
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

sunindextype SUNMatrixWrapper::num_nonzeros() const {
    if (!matrix_)
        return 0;
    
    if (SUNMatGetID(matrix_) == SUNMATRIX_SPARSE)
        return SM_INDEXPTRS_S(matrix_)[SM_NP_S(matrix_)];
    else
        throw std::domain_error("Non-zeros property only available for "
                                "sparse matrices");
}

sunindextype *SUNMatrixWrapper::indexvals() const {
    return indexvals_ptr_;
}

sunindextype *SUNMatrixWrapper::indexptrs() const {
    return indexptrs_ptr_;
}

int SUNMatrixWrapper::sparsetype() const {
    if (!matrix_)
        throw std::runtime_error("Cannot determine type of uninitialized "
                                 "matrices");
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
        int nonzeros_ = static_cast<int>(num_nonzeros());
        for (int i = 0; i < nonzeros_; ++i)
            data_ptr_[i] *= a;
    }
}

void SUNMatrixWrapper::multiply(N_Vector c, const_N_Vector b) const {
    multiply(gsl::make_span<realtype>(NV_DATA_S(c), NV_LENGTH_S(c)),
             gsl::make_span<const realtype>(NV_DATA_S(b), NV_LENGTH_S(b)));
}

static void check_dim(sunindextype n, sunindextype m, const std::string &name_n,
                      const std::string &name_m, const std::string &name_mat_n,
                      const std::string &name_mat_m) {
    if (n != m)
        throw std::invalid_argument("Dimension mismatch between number of "
                                    + name_n + " in " + name_mat_n + " ("
                                    + std::to_string(n)
                                    + ") and number of "
                                    + name_m + " in " + name_mat_m + " ("
                                    + std::to_string(m) + ")");
}

static void check_csc(const SUNMatrixWrapper *mat,
                      const std::string &fun,
                      const std::string &name_mat) {
    if (mat->matrix_id() != SUNMATRIX_SPARSE)
        throw std::invalid_argument(fun + " only implemented for "
                                    "sparse matrices, but "
                                    + name_mat + " is not sparse.");

    if (mat->sparsetype() != CSC_MAT)
        throw std::invalid_argument(fun + " only implemented for "
                                    "matrix type CSC, but "
                                    + name_mat + "is not of type CSC.");
}

void SUNMatrixWrapper::multiply(gsl::span<realtype> c,
                                gsl::span<const realtype> b) const {
    if (!matrix_)
        return;

    sunindextype nrows = rows();
    sunindextype ncols = columns();

    check_dim(nrows, c.size(), "rows", "elements", "A", "c");
    check_dim(ncols, b.size(), "cols", "elements", "A", "b");

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
        check_dim(cols.size() , c.size(), "selected columns",
                  "elements", "A", "c");
        check_dim(nrows, b.size(), "rows", "elements", "A", "b");
    } else {
        check_dim(nrows, c.size(), "rows", "elements", "A", "c");
        check_dim(ncols, b.size(), "columns", "elements", "A", "b");
    }

    check_csc(this, "Reordered multiply", "A");
    
    if (!num_nonzeros())
        return;

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


void SUNMatrixWrapper::sparse_multiply(SUNMatrixWrapper &C,
                                       const SUNMatrixWrapper &B) const {
    if (!matrix_ || !B.matrix_ || !C.matrix_)
        return;

    sunindextype nrows = rows();
    sunindextype ncols = columns();

    check_csc(this, "sparse_multiply", "A");
    check_csc(&B, "sparse_multiply", "B");
    check_csc(&C, "sparse_multiply", "C");

    check_dim(nrows, C.rows(), "rows", "rows", "A", "C");
    check_dim(C.columns(), B.columns(), "columns", "columns", "C", "B");
    check_dim(B.rows(), ncols, "rows", "columns", "B", "A");
    
    if (ncols == 0 || nrows == 0 || B.columns() == 0)
        return; // matrix will also have zero size
    
    if (num_nonzeros() == 0 || B.num_nonzeros() == 0)
        return; // nothing to multiply
    

    /* see https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/master/CSparse/Source/cs_multiply.c
     * modified such that we don't need to use CSparse memory structure and can
     * work with preallocated C. This should minimize number of necessary
     * reallocations as we can assume that C doesn't change size.
     */
    
    sunindextype nnz = 0; // this keeps track of the nonzero index in C
    auto Bx = B.data();
    auto Bi = B.indexvals();
    auto Bp = B.indexptrs();
    
    sunindextype bCol;
    sunindextype bcolptr;
    sunindextype ccolptr;
    auto Cx = C.data();
    auto Ci = C.indexvals();
    auto Cp = C.indexptrs();
    
    sunindextype m = nrows;
    sunindextype n = B.columns();
    
    auto w = std::vector<sunindextype>(m); // sparsity of C(:,j)
    auto x = std::vector<realtype>(m); // entries in C(:,j)
    
    for (bCol = 0; bCol < n; bCol++) // k in C(i,j) = sum_k A(i,k)*B(k,j)
    {
        Cp[bCol] = nnz;              /* column j of C starts here */
        if ((Bp[bCol+1] > Bp[bCol]) && (nnz + m > C.capacity()))
        {
            /*
             * if memory usage becomes a concern, remove the factor two here,
             * as it effectively trades memory efficiency against less
             * reallocations
             */
            C.reallocate(2*C.capacity() + m);
            // all pointers will change after reallocation
            Cx = C.data();
            Ci = C.indexvals();
            Cp = C.indexptrs();
        }
        for (bcolptr = Bp[bCol]; bcolptr < Bp[bCol+1]; bcolptr++)
        {
            nnz = scatter(Bi[bcolptr], Bx[bcolptr], w.data(), gsl::make_span(x),
                          bCol+1, &C, nnz);
            assert(nnz - Cp[bCol] <= m);
        }
        for (ccolptr = Cp[bCol]; ccolptr < nnz; ccolptr++)
            Cx[ccolptr] = x[Ci[ccolptr]]; // copy data to C
    }
    Cp[n] = nnz;
    assert(nnz <= C.capacity());
    /*
     * do not reallocate here since we rather keep a matrix that is a bit
     * bigger than repeatedly resizing this matrix.
     */
}

void SUNMatrixWrapper::sparse_add(const SUNMatrixWrapper &A, realtype alpha,
                                  const SUNMatrixWrapper &B, realtype beta) {
    // matrix_ == nullptr is allowed on the first call
    if (!A.matrix_ || !B.matrix_)
        return;

    sunindextype nrows = rows();
    sunindextype ncols = columns();

    check_csc(this, "sparse_multiply", "C");
    check_csc(&A, "sparse_multiply", "A");
    check_csc(&B, "sparse_multiply", "B");

    check_dim(nrows, A.rows(), "rows", "rows", "C", "A");
    check_dim(nrows, B.rows(), "rows", "rows", "C", "B");
    check_dim(ncols, A.columns(), "columns", "columns", "C", "A");
    check_dim(ncols, B.columns(), "columns", "columns", "C", "B");
    
    if (ncols == 0 || nrows == 0 ||
        (A.num_nonzeros() + B.num_nonzeros() == 0))
        return; // nothing to do
    

    /* see https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/master/CSparse/Source/cs_add.c
     * modified such that we don't need to use CSparse memory structure and can
     * work with preallocated C. This should minimize number of necessary
     * reallocations as we can assume that C doesn't change size.
     */
    
    sunindextype nnz = 0; // this keeps track of the nonzero index in C
    
    sunindextype cCol;
    sunindextype ccolptr;

    // first call, make sure that matrix is initialized with no capacity
    if(!capacity())
        reallocate(A.num_nonzeros() + B.num_nonzeros());
    auto Cx = data();
    auto Ci = indexvals();
    auto Cp = indexptrs();
    
    auto w = std::vector<sunindextype>(nrows);
    auto x = std::vector<realtype>(nrows);
    
    for (cCol = 0; cCol < ncols; cCol++)
    {
        Cp[cCol] = nnz;                          /* column j of C starts here */
        nnz = A.scatter(cCol, alpha, w.data(), gsl::make_span(x), cCol+1, this, nnz);
        nnz = B.scatter(cCol, beta, w.data(), gsl::make_span(x), cCol+1, this, nnz);
        // no reallocation should happen here
        for (ccolptr = Cp[cCol]; ccolptr < nnz; ccolptr++)
            Cx[ccolptr] = x[Ci[ccolptr]]; // copy data to C
    }
    Cp[ncols] = nnz;
    assert(nnz <= capacity());
    if (capacity() == A.num_nonzeros() + B.num_nonzeros())
        realloc(); // resize if necessary, will have correct size in future calls
}

sunindextype SUNMatrixWrapper::scatter(const sunindextype acol,
                                       const realtype beta,
                                       sunindextype *w,
                                       gsl::span<realtype> x,
                                       const sunindextype mark,
                                       SUNMatrixWrapper *C,
                                       sunindextype nnz) const {
    if (!matrix_)
        return nnz;
        
    check_csc(this, "scatter", "A");
    if (C && C->matrix_)
        check_csc(C, "scatter", "C");
    
    if (!num_nonzeros())
        return nnz;
    
    /* see https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/master/CSparse/Source/cs_scatter.c */
    
    sunindextype acolptr;
    sunindextype *Ci;
    if (C)
        Ci = C->indexvals();
    else
        Ci = nullptr;

    auto Ap = indexptrs();
    auto Ai = indexvals();
    auto Ax = data();
    for (acolptr = Ap[acol]; acolptr < Ap[acol+1]; acolptr++)
    {
        auto arow = Ai[acolptr];          /* A(arow,acol) is nonzero */
        assert(arow < static_cast<sunindextype>(x.size()));
        if (w && w[arow] < mark) {
            w[arow] = mark;               /* arow is new entry in C(:,*) */
            if (Ci)
                Ci[nnz++] = arow;         /* add arow to pattern of C(:,*) */
            x[arow] = beta * Ax[acolptr]; /* x(arow) = beta*A(arow,acol) */
        } else
            x[arow] += beta * Ax[acolptr];/* arow exists in C(:,*) already */
    }
    assert(!C || nnz <= C->capacity());
    return nnz;
}

// https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/master/CSparse/Source/cs_cumsum.c
/* p [0..n] = cumulative sum of c[0..n-1], and then copy p [0..n-1] into c */
static void cumsum(gsl::span<sunindextype> p, std::vector<sunindextype> &c) {
    sunindextype i;
    sunindextype nz = 0;
    assert(p.size() == c.size() + 1);
    for (i = 0; i < static_cast<sunindextype>(c.size()); i++)
    {
        p[i] = nz;
        nz += c[i];
        c[i] = p[i];             /* also copy p[0..n-1] back into c[0..n-1]*/
    }
    p[c.size()] = nz;
}

void SUNMatrixWrapper::transpose(SUNMatrix C, const realtype alpha,
                                 sunindextype blocksize) const{
    if (!matrix_ || !C)
        return;

    if (!((SUNMatGetID(C) == SUNMATRIX_SPARSE && SM_SPARSETYPE_S(C) == CSC_MAT)
          || SUNMatGetID(C) == SUNMATRIX_DENSE))
        throw std::domain_error("Not Implemented.");
    
    check_csc(this, "transpose", "A");
    if (SUNMatGetID(matrix_) == SUNMATRIX_SPARSE) {
        check_dim(rows(), SM_COLUMNS_S(C), "rows", "columns", "A", "C");
        check_dim(columns(), SM_ROWS_S(C), "columns", "rows", "A", "C");
        if (num_nonzeros() > SM_NNZ_S(C))
            std::invalid_argument("C must be allocated such that it can hold "
                                  "all nonzero values from A. Requires "
                                  + std::to_string(num_nonzeros()) + " was "
                                  + std::to_string(SM_NNZ_S(C)) + ".");
    } else {
        check_dim(rows(), SM_COLUMNS_D(C), "rows", "columns", "A", "C");
        check_dim(columns(), SM_ROWS_D(C), "columns", "rows", "A", "C");
    }
    
    auto ncols = columns();
    auto nrows = rows();
    
    assert(ncols % blocksize == 0);
    assert(nrows % blocksize == 0);
    
    if (!num_nonzeros() || !ncols || !nrows)
        return;
    
    // see https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/master/CSparse/Source/cs_transpose.c
    
    sunindextype aidx;
    sunindextype cidx;
    sunindextype acol;
    
    realtype *Cx;
    sunindextype *Ci;
    sunindextype *Cp;
    std::vector<sunindextype> w;
    
    sunindextype ccol;
    sunindextype crow;
    
    auto Ax = data();
    auto Ai = indexvals();
    auto Ap = indexptrs();
    
    if (SUNMatGetID(matrix_) == SUNMATRIX_SPARSE) {
        Cx = SM_DATA_S(C);
        Ci = SM_INDEXVALS_S(C);
        Cp = SM_INDEXPTRS_S(C);
        w = std::vector<sunindextype>(ncols);
        for (acol = 0; acol < nrows; acol++)                /* row counts */
            for (aidx = Ap[acol]; aidx < Ap[acol+1]; aidx++)
                w[(acol/blocksize)*blocksize + Ai[aidx] % blocksize]++;
        cumsum(gsl::make_span(Cp, ncols+1), w);             /* row pointers */
    }
    
    for (acol = 0; acol < nrows; acol++)
    {
        for (aidx = Ap[acol]; aidx < Ap[acol+1]; aidx++)
        {
            ccol = (acol/blocksize)*blocksize + Ai[aidx] % blocksize;
            crow = (Ai[aidx]/blocksize)*blocksize + acol % blocksize;
            if (SUNMatGetID(matrix_) == SUNMATRIX_SPARSE) {
                Ci[cidx = w[ccol]++] = crow;  /* place A(i,j) as entry C(j,i) */
                Cx[cidx] = alpha * Ax[aidx];
            } else {
                SM_ELEMENT_D(C, crow, ccol) = alpha * Ax[aidx];
            }
        }
    }
}

void SUNMatrixWrapper::to_dense(SUNMatrix D) const {
    if (!matrix_ || !D)
        return;
    check_csc(this, "to_dense", "A");
    check_dim(rows(), SM_ROWS_D(D), "rows", "rows", "A", "D");
    check_dim(columns(), SM_COLUMNS_D(D), "columns", "columns", "A", "D");
    
    SUNMatZero(D);
    if (!num_nonzeros())
        return;
        
    sunindextype icol;
    sunindextype idx;
    for (icol = 0; icol < columns(); ++icol)
        for (idx = indexptrs()[icol]; idx < indexptrs()[icol+1]; ++idx)
            SM_ELEMENT_D(D, indexvals()[idx], icol) = data()[idx];
}

void SUNMatrixWrapper::to_diag(N_Vector v) const {
    if (!matrix_ || !v)
        return;
    check_csc(this, "to_dense", "S");
    check_dim(rows(), columns(), "rows", "columns", "A", "A");
    check_dim(rows(), NV_LENGTH_S(v), "rows", "elements", "S", "v");
    
    N_VConst(0.0, v);
    if (!num_nonzeros())
        return;
        
    sunindextype icol;
    sunindextype idx;
    for (icol = 0; icol < columns(); ++icol)
        for (idx = indexptrs()[icol]; idx < indexptrs()[icol+1]; ++idx)
            if (indexvals()[idx] == icol)
                NV_Ith_S(v, icol) = data()[idx];
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

