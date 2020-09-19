#include <amici/sundials_matrix_wrapper.h>
#include <sundials/sundials_matrix.h> // return codes

#include <amici/cblas.h>

#include <new> // bad_alloc
#include <utility>
#include <stdexcept> // invalid_argument and domain_error
#include <assert.h>

namespace amici {

SUNMatrixWrapper::SUNMatrixWrapper(sunindextype M, sunindextype N,
                                   sunindextype NNZ, int sparsetype)
    : matrix_(SUNSparseMatrix(M, N, NNZ, sparsetype)),
    id_(SUNMatGetID(matrix_))  {

    if (sparsetype != CSC_MAT && sparsetype != CSR_MAT)
        throw std::invalid_argument("Invalid sparsetype. Must be CSC_MAT or "
                                    "CSR_MAT");

    if (NNZ && M && N && !matrix_)
        throw std::bad_alloc();
        
    assert(num_nonzeros() == 0);
    assert(NNZ == capacity());
    assert(M == rows() || !matrix_);
    assert(N == columns() || !matrix_);
}

SUNMatrixWrapper::SUNMatrixWrapper(sunindextype M, sunindextype N)
    : matrix_(SUNDenseMatrix(M, N)),
    id_(SUNMatGetID(matrix_)) {
    if (M && N && !matrix_)
        throw std::bad_alloc();
        
    assert(M == rows());
    assert(N == columns());
}

SUNMatrixWrapper::SUNMatrixWrapper(sunindextype M, sunindextype ubw,
                                   sunindextype lbw)
    : matrix_(SUNBandMatrix(M, ubw, lbw)),
    id_(SUNMatGetID(matrix_)) {
    if (M && !matrix_)
        throw std::bad_alloc();
}

SUNMatrixWrapper::SUNMatrixWrapper(const SUNMatrixWrapper &A, realtype droptol,
                                   int sparsetype)
    : id_(A.id_){
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
}

SUNMatrixWrapper::SUNMatrixWrapper(SUNMatrix mat) : matrix_(mat),
    id_(SUNMatGetID(mat)), ownmat(false) {}

SUNMatrixWrapper::~SUNMatrixWrapper() {
    if (matrix_ && ownmat)
        SUNMatDestroy(matrix_);
}

SUNMatrixWrapper::SUNMatrixWrapper(const SUNMatrixWrapper &other)
    : id_(other.id_) {
    if (!other.matrix_)
        return;

    matrix_ = SUNMatClone(other.matrix_);
    if (!matrix_)
        throw std::bad_alloc();

    SUNMatCopy(other.matrix_, matrix_);
}

SUNMatrixWrapper::SUNMatrixWrapper(SUNMatrixWrapper &&other) : id_(other.id_) {
    std::swap(matrix_, other.matrix_);
}

SUNMatrixWrapper &SUNMatrixWrapper::operator=(const SUNMatrixWrapper &other) {
    if(&other == this)
        return *this;
    return *this = SUNMatrixWrapper(other);
}

SUNMatrixWrapper &SUNMatrixWrapper::operator=(SUNMatrixWrapper &&other) {
    std::swap(matrix_, other.matrix_);
    id_ = other.id_;
    return *this;
}

void SUNMatrixWrapper::reallocate(sunindextype NNZ) {
    if (sparsetype() != CSC_MAT && sparsetype() != CSR_MAT)
        throw std::invalid_argument("Invalid sparsetype. Must be CSC_MAT or "
                                    "CSR_MAT.");
    
    if (int ret = SUNSparseMatrix_Reallocate(matrix_, NNZ) != SUNMAT_SUCCESS)
        throw std::runtime_error("SUNSparseMatrix_Reallocate failed with "
                                 "error code " + std::to_string(ret) + ".");

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
    
    assert(capacity() ^ !matrix_);
    assert(capacity() == num_nonzeros());
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

static void check_sparse(const SUNMatrix_ID id) {
    if (id != SUNMATRIX_SPARSE)
        throw std::domain_error("Only available for sparse matrices");
}

static void check_dense(const SUNMatrix_ID id) {
    if (id != SUNMATRIX_DENSE)
        throw std::domain_error("Only available for dense matrices");
}


realtype SUNMatrixWrapper::get_data(sunindextype idx) const{
    assert(matrix_);
    check_sparse(matrix_id());
    assert(idx < SM_NNZ_S(matrix_));
    return SM_DATA_S(matrix_)[idx];
};

realtype SUNMatrixWrapper::get_data(sunindextype irow, sunindextype icol) const{
    assert(matrix_);
    check_dense(matrix_id());
    assert(irow < rows());
    assert(icol < columns());
    return SM_ELEMENT_D(matrix_, irow, icol);
};


void SUNMatrixWrapper::set_data(sunindextype idx, realtype data) {
    assert(matrix_);
    check_sparse(matrix_id());
    assert(idx < SM_NNZ_S(matrix_));
    SM_DATA_S(matrix_)[idx] = data;
}

void SUNMatrixWrapper::set_data(sunindextype irow, sunindextype icol,
                                realtype data) {
    assert(matrix_);
    check_dense(matrix_id());
    assert(irow < rows());
    assert(icol < columns());
    SM_ELEMENT_D(matrix_, irow, icol) = data;
}

const realtype *SUNMatrixWrapper::data() const {
    if (!matrix_)
        return nullptr;
        
    switch (SUNMatGetID(matrix_)) {
    case SUNMATRIX_DENSE:
        return SM_DATA_D(matrix_);
    case SUNMATRIX_SPARSE:
        return SM_DATA_S(matrix_);
    default:
        throw std::domain_error("Not Implemented.");
    }
}

realtype *SUNMatrixWrapper::data() {
    return const_cast<realtype *>(
        const_cast<const SUNMatrixWrapper *> (this)->data());
}

sunindextype SUNMatrixWrapper::get_indexval(sunindextype idx) const {
    assert(matrix_);
    check_sparse(matrix_id());
    assert(idx < SM_NNZ_S(matrix_));
    return SM_INDEXVALS_S(matrix_)[idx];
}

void SUNMatrixWrapper::set_indexval(sunindextype idx, sunindextype val) {
    assert(matrix_);
    check_sparse(matrix_id());
    assert(idx < SM_NNZ_S(matrix_));
    SM_INDEXVALS_S(matrix_)[idx] = val;
}

void SUNMatrixWrapper::set_indexvals(const gsl::span<const sunindextype> vals) {
    assert(matrix_);
    check_sparse(matrix_id());
    assert(static_cast<sunindextype>(vals.size()) == SM_NNZ_S(matrix_));
    std::copy_n(vals.begin(), SM_NNZ_S(matrix_), SM_INDEXVALS_S(matrix_));
}

sunindextype SUNMatrixWrapper::get_indexptr(sunindextype ptr_idx) const {
    assert(matrix_);
    check_sparse(matrix_id());
    assert(ptr_idx <= SM_NP_S(matrix_));
    return SM_INDEXPTRS_S(matrix_)[ptr_idx];
}

void SUNMatrixWrapper::set_indexptr(sunindextype ptr_idx, sunindextype ptr) {
    assert(matrix_);
    check_sparse(matrix_id());
    assert(ptr_idx <= SM_NP_S(matrix_));
    assert(ptr <= capacity());
    SM_INDEXPTRS_S(matrix_)[ptr_idx] = ptr;
}

void SUNMatrixWrapper::set_indexptrs(const gsl::span<const sunindextype> ptrs) {
    assert(matrix_);
    check_sparse(matrix_id());
    assert(static_cast<sunindextype>(ptrs.size()) == SM_NP_S(matrix_) + 1);
    std::copy_n(ptrs.begin(), SM_NP_S(matrix_) + 1, SM_INDEXPTRS_S(matrix_));
}

int SUNMatrixWrapper::sparsetype() const {
    if (!matrix_)
        throw std::runtime_error("Cannot determine type of uninitialized "
                                 "matrices");
    if (SUNMatGetID(matrix_) == SUNMATRIX_SPARSE)
        return SM_SPARSETYPE_S(matrix_);
    throw std::domain_error("Function only available for sparse matrices");
}

void SUNMatrixWrapper::scale(realtype a) {
    if (matrix_) {
        auto nonzeros_ = SM_NNZ_S(matrix_);
        for (sunindextype i = 0; i < nonzeros_; ++i)
            SM_DATA_S(matrix_)[i] *= a;
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
            return;
        }
        check_csc(this, "multiply", "A");
        for (sunindextype icol = 0; icol < ncols; ++icol) {
            scatter(icol, b.at(icol), nullptr, c, icol+1, nullptr, 0);
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
    sunindextype idx;
    if (transpose) {
        for (int icols = 0; icols < (int)cols.size(); ++icols)
            for (idx = get_indexptr(cols.at(icols));
                 idx < get_indexptr(cols.at(icols) + 1); ++idx)
                c.at(icols) += get_data(idx) * b.at(get_indexval(idx));
    } else {
        for (sunindextype icols = 0; icols < ncols; ++icols)
            for (idx = get_indexptr(cols.at(icols));
                 idx < get_indexptr(cols.at(icols)+1); ++idx)
                c.at(get_indexval(idx)) += get_data(idx) * b.at(icols);
    }
}


void SUNMatrixWrapper::sparse_multiply(SUNMatrixWrapper &C,
                                       const SUNMatrixWrapper &B) const {
    if (!matrix_ || !B.matrix_ || !C.matrix_)
        return;

    sunindextype nrows = rows();
    sunindextype ncols = columns();
    sunindextype bcols = B.columns();

    check_csc(this, "sparse_multiply", "A");
    check_csc(&B, "sparse_multiply", "B");
    check_csc(&C, "sparse_multiply", "C");

    check_dim(nrows, C.rows(), "rows", "rows", "A", "C");
    check_dim(C.columns(), B.columns(), "columns", "columns", "C", "B");
    check_dim(B.rows(), ncols, "rows", "columns", "B", "A");
    
    C.zero();
    
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
    sunindextype bcol;
    sunindextype bidx;
    sunindextype cidx;
    
    auto w = std::vector<sunindextype>(nrows); // sparsity of C(:,j)
    auto x = std::vector<realtype>(nrows); // entries in C(:,j)
    
    for (bcol = 0; bcol < bcols; bcol++) // k in C(i,j) = sum_k A(i,k)*B(k,j)
    {
        C.set_indexptr(bcol, nnz);              /* column j of C starts here */
        if ((B.get_indexptr(bcol+1) > B.get_indexptr(bcol))
            && (nnz + nrows > C.capacity()))
        {
            /*
             * if memory usage becomes a concern, remove the factor two here,
             * as it effectively trades memory efficiency against less
             * reallocations
             */
            C.reallocate(2*C.capacity() + nrows);
        }
        for (bidx = B.get_indexptr(bcol); bidx < B.get_indexptr(bcol+1); bidx++)
        {
            nnz = scatter(B.get_indexval(bidx), B.get_data(bidx),
                          w.data(), gsl::make_span(x), bcol+1, &C, nnz);
            assert(nnz - C.get_indexptr(bcol) <= nrows);
        }
        for (cidx = C.get_indexptr(bcol); cidx < nnz; cidx++)
            C.set_data(cidx, x.at(C.get_indexval(cidx))); // copy data to C
    }
    C.set_indexptr(C.columns(), nnz);
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
    
    zero();
    
    if (ncols == 0 || nrows == 0 ||
        (A.num_nonzeros() + B.num_nonzeros() == 0))
        return; // nothing to do
    

    /* see https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/master/CSparse/Source/cs_add.c
     * modified such that we don't need to use CSparse memory structure and can
     * work with preallocated C. This should minimize number of necessary
     * reallocations as we can assume that C doesn't change size.
     */
    
    sunindextype nnz = 0; // this keeps track of the nonzero index in C
    
    sunindextype ccol;
    sunindextype cidx;

    // first call, make sure that matrix is initialized with no capacity
    if(!capacity())
        reallocate(A.num_nonzeros() + B.num_nonzeros());

    auto w = std::vector<sunindextype>(nrows);
    auto x = std::vector<realtype>(nrows);
    
    for (ccol = 0; ccol < ncols; ccol++)
    {
        set_indexptr(ccol, nnz);                          /* column j of C starts here */
        nnz = A.scatter(ccol, alpha, w.data(), gsl::make_span(x), ccol+1, this,
                        nnz);
        nnz = B.scatter(ccol, beta, w.data(), gsl::make_span(x), ccol+1, this,
                        nnz);
        // no reallocation should happen here
        for (cidx = get_indexptr(ccol); cidx < nnz; cidx++) {
            set_data(cidx, x.at(get_indexval(cidx))); // copy data to C
        }
    }
    set_indexptr(ncols, nnz);
    if (capacity() == A.num_nonzeros() + B.num_nonzeros())
        realloc(); // resize if necessary, will have correct size in future calls
}

void SUNMatrixWrapper::sparse_sum(const std::vector<SUNMatrixWrapper> &mats) {
    // matrix_ == nullptr is allowed on the first call
    if (std::all_of(mats.begin(), mats.end(), [](const SUNMatrixWrapper &m){return !m.matrix_;}))
        return;

    sunindextype nrows = rows();
    sunindextype ncols = columns();

    check_csc(this, "sparse_multiply", "A");
    int max_total_nonzero = 0;
    for (auto & mat : mats) {
        check_csc(&mat, "sparse_multiply", "mat");
        check_dim(nrows, mat.rows(), "rows", "rows", "A", "mat");
        check_dim(ncols, mat.columns(), "columns", "columns", "A", "mat");
        max_total_nonzero += mat.num_nonzeros();
    }
    
    zero();
    
    if (ncols == 0 || nrows == 0 || max_total_nonzero == 0)
        return; // nothing to do

    /* see https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/master/CSparse/Source/cs_add.c
     * modified such that we don't need to use CSparse memory structure and can
     * work with preallocated C. This should minimize number of necessary
     * reallocations as we can assume that C doesn't change size.
     */
    
    sunindextype nnz = 0; // this keeps track of the nonzero index in C
    
    sunindextype acol;
    sunindextype aidx;
    // first call, make sure that matrix is initialized with no capacity
    if(!capacity())
        reallocate(max_total_nonzero);
    
    auto w = std::vector<sunindextype>(nrows);
    auto x = std::vector<realtype>(nrows);
    
    for (acol = 0; acol < ncols; acol++)
    {
        set_indexptr(acol, nnz);                       /* column j of A starts here */
        for (auto & mat : mats)
            nnz = mat.scatter(acol, 1.0, w.data(), gsl::make_span(x), acol+1,
                              this, nnz);
        // no reallocation should happen here
        for (aidx = get_indexptr(acol); aidx < nnz; aidx++) {
            set_data(aidx, x.at(get_indexval(aidx))); // copy data to C
        }
    }
    set_indexptr(ncols, nnz);
    if (capacity() == max_total_nonzero)
        realloc(); // resize if necessary
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
    
    sunindextype aidx;
    for (aidx = get_indexptr(acol); aidx < get_indexptr(acol+1); aidx++)
    {
        auto arow = get_indexval(aidx);          /* A(arow,acol) is nonzero */
        if (w && w[arow] < mark) {
            w[arow] = mark;                      /* arow is new entry in C(:,*) */
            if (C)
                C->set_indexval(nnz++, arow);    /* add arow to pattern of C(:,*) */
            x.at(arow) = beta * get_data(aidx);  /* x(arow) = beta*A(arow,acol) */
        } else
            x.at(arow) += beta * get_data(aidx); /* arow exists in C(:,*) already */
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

void SUNMatrixWrapper::transpose(SUNMatrixWrapper &C, const realtype alpha,
                                 sunindextype blocksize) const{
    if (!matrix_ || !C.matrix_)
        return;

    if (!((C.matrix_id() == SUNMATRIX_SPARSE && C.sparsetype() == CSC_MAT)
          || C.matrix_id() == SUNMATRIX_DENSE))
        throw std::domain_error("Not Implemented.");
    
    check_csc(this, "transpose", "A");
    check_dim(rows(), C.rows(), "rows", "columns", "A", "C");
    check_dim(columns(), C.columns(), "columns", "rows", "A", "C");
    if (C.matrix_id() == SUNMATRIX_SPARSE) {
        if (!C.capacity() && num_nonzeros())
            C.reallocate(num_nonzeros());
        if (num_nonzeros() > C.capacity())
            std::invalid_argument("C must be allocated such that it can hold "
                                  "all nonzero values from A. Requires "
                                  + std::to_string(num_nonzeros()) + " was "
                                  + std::to_string(C.capacity()) + ".");
    }

    auto ncols = columns();
    auto nrows = rows();
    
    assert(ncols % blocksize == 0);
    assert(nrows % blocksize == 0);
    C.zero();
    
    if (!num_nonzeros() || !ncols || !nrows)
        return;
    
    // see https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/master/CSparse/Source/cs_transpose.c
    
    sunindextype aidx;
    sunindextype cidx;
    sunindextype acol;

    std::vector<sunindextype> w;
    
    sunindextype ccol;
    sunindextype crow;
    sunindextype widx;
    
    if (C.matrix_id()== SUNMATRIX_SPARSE) {
        w = std::vector<sunindextype>(ncols);
        for (acol = 0; acol < nrows; acol++)                /* row counts */
            for (aidx = get_indexptr(acol);
                 aidx < get_indexptr(acol+1); aidx++) {
                widx = (acol/blocksize)*blocksize + get_indexval(aidx) % blocksize;
                w.at(widx)++;
                assert(w[widx] <= nrows);
            }
        /* row pointers */
        cumsum(gsl::make_span(SM_INDEXPTRS_S(C.matrix_), C.columns()+1), w);
    }
    
    for (acol = 0; acol < nrows; acol++)
    {
        for (aidx = get_indexptr(acol); aidx < get_indexptr(acol+1); aidx++)
        {
            ccol = (acol/blocksize)*blocksize + get_indexval(aidx) % blocksize;
            crow = (get_indexval(aidx)/blocksize)*blocksize + acol % blocksize;
            assert(crow < nrows);
            assert(ccol < ncols);
            if (C.matrix_id() == SUNMATRIX_SPARSE) {
                assert(aidx < capacity());
                cidx = w.at(ccol)++;
                C.set_indexval(cidx, crow);  /* place A(i,j) as entry C(j,i) */
                C.set_data(cidx, alpha * get_data(aidx));
            } else {
                C.set_data(crow, ccol, alpha * get_data(aidx));
            }
        }
    }
}

void SUNMatrixWrapper::to_dense(SUNMatrixWrapper &D) const {
    if (!matrix_ || !D.matrix_)
        return;
    check_csc(this, "to_dense", "A");
    check_dim(rows(), D.rows(), "rows", "rows", "A", "D");
    check_dim(columns(), D.columns(), "columns", "columns", "A", "D");
    
    D.zero();
    if (!num_nonzeros())
        return;
        
    sunindextype icol;
    sunindextype idx;
    for (icol = 0; icol < columns(); ++icol)
        for (idx = get_indexptr(icol); idx < get_indexptr(icol+1); ++idx) {
            D.set_data(get_indexval(idx), icol, get_data(idx));
        }
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
        for (idx = get_indexptr(icol); idx < get_indexptr(icol+1); ++idx)
            if (get_indexval(idx) == icol)
                NV_Ith_S(v, icol) = get_data(idx);
}


void SUNMatrixWrapper::zero()
{
    if (!matrix_)
        return;
    if(int res = SUNMatZero(matrix_))
        throw std::runtime_error("SUNMatrixWrapper::zero() failed with "
                                 + std::to_string(res) + ".");
}

const SUNMatrix SUNMatrixWrapper::get() const { return matrix_; }

} // namespace amici

