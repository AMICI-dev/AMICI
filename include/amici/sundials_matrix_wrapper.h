#ifndef AMICI_SUNDIALS_MATRIX_WRAPPER_H
#define AMICI_SUNDIALS_MATRIX_WRAPPER_H

#include <sunmatrix/sunmatrix_band.h>   // SUNMatrix_Band
#include <sunmatrix/sunmatrix_dense.h>  // SUNMatrix_Dense
#include <sunmatrix/sunmatrix_sparse.h> // SUNMatrix_Sparse

#include <gsl/gsl-lite.hpp>

#include <vector>

#include "amici/vector.h"

namespace amici {

/**
 * @brief A RAII wrapper for SUNMatrix structs.
 *
 * This can create dense, sparse, or banded matrices using the respective
 * constructor.
 */
class SUNMatrixWrapper {
  public:
    SUNMatrixWrapper() = default;

    /**
     * @brief Create sparse matrix. See SUNSparseMatrix in sunmatrix_sparse.h
     * @param M Number of rows
     * @param N Number of columns
     * @param NNZ Number of nonzeros
     * @param sparsetype Sparse type
     */
    SUNMatrixWrapper(sunindextype M, sunindextype N, sunindextype NNZ,
                     int sparsetype);

    /**
     * @brief Create dense matrix. See SUNDenseMatrix in sunmatrix_dense.h
     * @param M Number of rows
     * @param N Number of columns
     */
    SUNMatrixWrapper(sunindextype M, sunindextype N);

    /**
     * @brief Create banded matrix. See SUNBandMatrix in sunmatrix_band.h
     * @param M Number of rows and columns
     * @param ubw Upper bandwidth
     * @param lbw Lower bandwidth
     */
    SUNMatrixWrapper(sunindextype M, sunindextype ubw, sunindextype lbw);

    /**
     * @brief Create sparse matrix from dense or banded matrix. See
     * SUNSparseFromDenseMatrix and SUNSparseFromBandMatrix in
     * sunmatrix_sparse.h
     * @param A Wrapper for dense matrix
     * @param droptol tolerance for dropping entries
     * @param sparsetype Sparse type
     */
    SUNMatrixWrapper(const SUNMatrixWrapper &A, realtype droptol,
                     int sparsetype);

    /**
     * @brief Wrap existing SUNMatrix
     * @param mat
     */
    explicit SUNMatrixWrapper(SUNMatrix mat);

    ~SUNMatrixWrapper();

    /**
     * @brief Copy constructor
     * @param other
     */
    SUNMatrixWrapper(const SUNMatrixWrapper &other);

    /**
     * @brief Move constructor
     * @param other
     */
    SUNMatrixWrapper(SUNMatrixWrapper &&other);

    /**
     * @brief Copy assignment
     * @param other
     * @return
     */
    SUNMatrixWrapper &operator=(const SUNMatrixWrapper &other);

    /**
     * @brief Move assignment
     * @param other
     * @return
     */
    SUNMatrixWrapper &operator=(SUNMatrixWrapper &&other);
    
    /**
     * @brief Reallocate space for sparse matrix according to specified nnz
     * @param nnz new number of nonzero entries
     */
    void reallocate(sunindextype nnz);
    
    /**
     * @brief Reallocate space for sparse matrix to used space according to last entry in indexptrs
     */
    void realloc();

    /**
     * @brief Get the wrapped SUNMatrix
     * @return raw SunMatrix object
     * @note Even though the returned matrix_ pointer is const qualified, matrix_->content will not be const.
     * This is a shortcoming in the underlying C library, which we cannot address and it is not intended that
     * any of those values are modified externally. If matrix_->content is manipulated,
     * cpp:meth:SUNMatrixWrapper:`refresh` needs to be called.
     */
    const SUNMatrix get() const;

    /**
     * @brief Get the number of rows
     * @return number of rows
     */
    sunindextype rows() const;

    /**
     * @brief Get the number of columns
     * @return number of columns
     */
    sunindextype columns() const;

    /**
     * @brief Get the number of specified non-zero elements (sparse matrices only)
     * @note value will be 0 before indexptrs are set.
     * @return number of nonzero entries
     */
    sunindextype num_nonzeros() const;
    
    /**
     * @brief Get the number of indexptrs that can be specified (sparse matrices only)
     * @return number of indexptrs
     */
    sunindextype num_indexptrs() const;
    
    /**
     * @brief Get the number of allocated data elements
     * @return number of allocated entries
     */
    sunindextype capacity() const;
    
    /**
     * @brief Get  raw data of a sparse matrix
     * @return pointer to first data entry
     */
    realtype *data();
    
    /**
     * @brief Get const raw data of a sparse matrix
     * @return pointer to first data entry
     */
    const realtype *data() const;
    
    /**
     * @brief Get data of a sparse matrix
     * @param idx data index
     * @return idx-th data entry
     */
    realtype get_data(sunindextype idx) const;
    
    /**
     * @brief Get data entry for a dense matrix
     * @param irow row
     * @param icol col
     * @return A(irow,icol)
     */
    realtype get_data(sunindextype irow, sunindextype icol) const;
    
    /**
     * @brief Set data entry for a sparse matrix
     * @param idx data index
     * @param data data for idx-th entry
     */
    void set_data(sunindextype idx, realtype data);
    
    /**
     * @brief Set data entry for a dense matrix
     * @param irow row
     * @param icol col
     * @param data data for idx-th entry
     */
    void set_data(sunindextype irow, sunindextype icol, realtype data);

    /**
     * @brief Get the index value of a sparse matrix
     * @param idx data index
     * @return row (CSC) or column (CSR) for idx-th data entry
     */
    sunindextype get_indexval(sunindextype idx) const;
    
    /**
     * @brief Set the index value of a sparse matrix
     * @param idx data index
     * @param val row (CSC) or column (CSR) for idx-th data entry
     */
    void set_indexval(sunindextype idx, sunindextype val);
    
    /**
     * @brief Set the index values of a sparse matrix
     * @param vals rows (CSC) or columns (CSR) for data entries
     */
    void set_indexvals(const gsl::span<const sunindextype> vals);
    
    /**
     * @brief Get the index pointer of a sparse matrix
     * @param ptr_idx pointer index
     * @return index where the ptr_idx-th column (CSC) or row (CSR) starts
     */
    sunindextype get_indexptr(sunindextype ptr_idx) const;
    
    /**
     * @brief Set the index pointer of a sparse matrix
     * @param ptr_idx pointer index
     * @param ptr data-index where the ptr_idx-th column (CSC) or row (CSR) starts
     */
    void set_indexptr(sunindextype ptr_idx, sunindextype ptr);
    
    /**
     * @brief Set the index pointers of a sparse matrix
     * @param ptrs starting data-indices where the columns (CSC) or rows (CSR) start
     */
    void set_indexptrs(const gsl::span<const sunindextype> ptrs);

    /**
     * @brief Get the type of sparse matrix
     * @return matrix type
     */
    int sparsetype() const;

    /**
     * @brief multiply with a scalar (in-place)
     * @param a scalar value to multiply matrix
     */
    void scale(realtype a);

    /**
     * @brief N_Vector interface for multiply
     * @param c output vector, may already contain values
     * @param b multiplication vector
     */
    void multiply(N_Vector c, const_N_Vector b) const;

    /**
     * @brief Perform matrix vector multiplication c += A*b
     * @param c output vector, may already contain values
     * @param b multiplication vector
     */
    void multiply(gsl::span<realtype> c, gsl::span<const realtype> b) const;

    /**
     * @brief Perform reordered matrix vector multiplication c += A[:,cols]*b
     * @param c output vector, may already contain values
     * @param b multiplication vector
     * @param cols int vector for column reordering
     * @param transpose bool transpose A before multiplication
     */
    void multiply(N_Vector c,
                  const N_Vector b,
                  gsl::span <const int> cols,
                  bool transpose) const;

    /**
     * @brief Perform reordered matrix vector multiplication c += A[:,cols]*b
     * @param c output vector, may already contain values
     * @param b multiplication vector
     * @param cols int vector for column reordering
     * @param transpose bool transpose A before multiplication
     */
    void multiply(gsl::span<realtype> c,
                  gsl::span<const realtype> b,
                  gsl::span <const int> cols,
                  bool transpose) const;

    /**
     * @brief Perform matrix matrix multiplication C = A * B for sparse A, B, C
     * @param C output matrix,
     * @param B multiplication matrix
     * @note will overwrite existing data, indexptrs, indexvals for C, but will use preallocated space for these vars
     */
    void sparse_multiply(SUNMatrixWrapper &C,
                         const SUNMatrixWrapper &B) const;
    
    /**
     * @brief Perform sparse matrix matrix addition C = alpha * A +  beta * B
     * @param A addition matrix
     * @param alpha scalar A
     * @param B addition matrix
     * @param beta scalar B
     * @note will overwrite existing data, indexptrs, indexvals for C, but will use preallocated space for these vars
     */
    void sparse_add(const SUNMatrixWrapper &A, realtype alpha,
                    const SUNMatrixWrapper &B, realtype beta);
    
    /**
     * @brief Perform matrix-matrix addition A = sum(mats(0)...mats(len(mats)))
     * @param mats vector of sparse matrices
     * @note will overwrite existing data, indexptrs, indexvals for A, but will use preallocated space for these vars
     */
    void sparse_sum(const std::vector<SUNMatrixWrapper> &mats);
    
    /**
     * @brief Compute x = x + beta * A(:,k), where x is a dense vector and A(:,k) is sparse, and update
     * the sparsity pattern for C(:,j) if applicable
     *
     * This function currently has two purposes:
     *   - perform parts of sparse matrix-matrix multiplication C(:,j)=A(:,k)*B(k,j)
     *    enabled by passing beta=B(k,j), x=C(:,j), C=C, w=sparsity of C(:,j) from B(k,0...j-1), nnz=nnz(C(:,0...j-1)
     *   - add the k-th column of the sparse matrix A multiplied by beta to the dense vector x.
     *    enabled by passing beta=*, x=x, C=nullptr, w=nullptr, nnz=*
     *
     * @param k column index
     * @param beta scaling factor
     * @param w index workspace, (w[i]<mark) indicates non-zeroness of C(i,j) (dimension: m),
     * if this is a nullptr, sparsity pattern of C will not be updated (if applicable).
     * @param x dense output vector (dimension: m)
     * @param mark marker for w to indicate nonzero pattern
     * @param C sparse output matrix, if this is a nullptr, sparsity pattern of C will not be updated
     * @param nnz number of nonzeros that were already written to C
     * @return updated number of nonzeros in C
     */
    sunindextype scatter(const sunindextype k, const realtype beta,
                         sunindextype *w, gsl::span<realtype> x,
                         const sunindextype mark,
                         SUNMatrixWrapper *C, sunindextype nnz) const;
    
    /**
     * @brief Compute transpose A' of sparse matrix A and writes it to the matrix C = alpha * A'
     *
     * @param C output matrix (sparse or dense)
     * @param alpha scalar multiplier
     * @param blocksize blocksize for transposition. For full matrix transpose set to ncols/nrows
     */
    void transpose(SUNMatrixWrapper &C, const realtype alpha,
                   sunindextype blocksize) const;
    
    /**
     * @brief Writes a sparse matrix A to a dense matrix D.
     *
     * @param D dense output matrix
     */
    void to_dense(SUNMatrixWrapper &D) const;
    
    /**
     * @brief Writes the diagonal of sparse matrix A to a dense vector v.
     *
     * @param v dense outut vector
     */
    void to_diag(N_Vector v) const;

    /**
     * @brief Set to 0.0, for sparse matrices also resets indexptr/indexvals
     */
    void zero();
    
    /**
     * @brief Get matrix id
     * @return SUNMatrix_ID
     */
    SUNMatrix_ID matrix_id() const {return id_;};
    
    /**
     * @brief Update internal cache, needs to be called after external manipulation of matrix_->content
     */
    void refresh();
    
  private:

    /**
     * @brief SUNMatrix to which all methods are applied
     */
    SUNMatrix matrix_ {nullptr};
    
    /**
     * @brief cache for SUNMatrixGetId(matrix_)
     */
    SUNMatrix_ID id_ {SUNMATRIX_CUSTOM};
    
    /**
     * @brief cache for SUNMatrixGetId(matrix_)
     */
    int sparsetype_ {CSC_MAT};
    
    /**
     * @brief cache for SM_INDEXPTRS_S(matrix_)[SM_NP_S(matrix_)]
     */
    sunindextype num_nonzeros_ {0};
    /**
     * @brief cache for SM_NNZ_S(matrix_)
     */
    sunindextype capacity_ {0};
    
    /**
     * @brief cache for SM_DATA_S(matrix_)
     */
    realtype *data_ {nullptr};
    /**
     * @brief cache for SM_INDEXPTRS_S(matrix_)
     */
    sunindextype *indexptrs_ {nullptr};
    /**
     * @brief cache for SM_INDEXVALS_S(matrix_)
     */
    sunindextype *indexvals_ {nullptr};
    
    /**
     * @brief cache for SM_ROWS_X(matrix_)
     */
    sunindextype num_rows_ {0};
    /**
     * @brief cache for SM_COLUMS_X(matrix_)
     */
    sunindextype num_columns_ {0};
    /**
     * @brief cache for SM_NP_S(matrix_)
     */
    sunindextype num_indexptrs_ {0};
    
    /**
     * @brief call update_ptrs & update_size
     */
    void finish_init();
    /**
     * @brief update data_, indexptrs_, indexvals_ if applicable
     */
    void update_ptrs();
    /**
     * @brief update num_rows_, num_columns_, num_indexptrs if applicable
     */
    void update_size();
    /**
     * @brief indicator whether this wrapper allocated matrix_ and is responsible for deallocation
     */
    bool ownmat = true;
};

} // namespace amici

namespace gsl {
/**
 * @brief Create span from SUNMatrix
 * @param m SUNMatrix
 * @return Created span
 */
inline span<realtype> make_span(SUNMatrix m)
{
    switch (SUNMatGetID(m)) {
    case SUNMATRIX_DENSE:
        return span<realtype>(SM_DATA_D(m), SM_LDATA_D(m));
    case SUNMATRIX_SPARSE:
        return span<realtype>(SM_DATA_S(m), SM_NNZ_S(m));
    default:
        throw amici::AmiException("Unimplemented SUNMatrix type for make_span");
    }
}
} // namespace gsl

#endif // AMICI_SUNDIALS_MATRIX_WRAPPER_H
