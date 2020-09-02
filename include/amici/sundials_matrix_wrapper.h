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
    SUNMatrixWrapper(int M, int N, int NNZ, int sparsetype);

    /**
     * @brief Create dense matrix. See SUNDenseMatrix in sunmatrix_dense.h
     * @param M Number of rows
     * @param N Number of columns
     */
    SUNMatrixWrapper(int M, int N);

    /**
     * @brief Create banded matrix. See SUNBandMatrix in sunmatrix_band.h
     * @param M Number of rows and columns
     * @param ubw Upper bandwidth
     * @param lbw Lower bandwidth
     */
    SUNMatrixWrapper(int M, int ubw, int lbw);

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
     * @brief Access raw data
     * @return raw data pointer
     */
    realtype *data() const;

    /**
     * @brief Get the wrapped SUNMatrix
     * @return SlsMat
     */
    SUNMatrix get() const;

    /**
     * @brief Get the number of rows
     * @return number
     */
    sunindextype rows() const;

    /**
     * @brief Get the number of columns
     * @return number
     */
    sunindextype columns() const;

    /**
     * @brief Get the number of specified non-zero elements (sparse matrices only)
     * @note value will be 0 before indexptrs are set.
     * @return number
     */
    sunindextype num_nonzeros() const;
    
    /**
     * @brief Get the number of allocated non-zero elements (sparse matrices only)
     * @return number
     */
    sunindextype capacity() const;

    /**
     * @brief Get the index values of a sparse matrix
     * @return index array
     */
    sunindextype *indexvals() const;

    /**
     * @brief Get the index pointers of a sparse matrix
     * @return index array
     */
    sunindextype *indexptrs() const;

    /**
     * @brief Get the type of sparse matrix
     * @return index array
     */
    int sparsetype() const;

    /**
     * @brief reset data to zeroes
     */
    void reset();

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
     * @brief Perform matrix matrix addition C = alpha * A +  beta * B
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
    void sparse_sum(const std::vector<SUNMatrixWrapper> mats);
    
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
     * @brief Set to 0.0
     */
    void zero();
    
    /**
     * @brief Get matrix id
     * @return SUNMatrix_ID
     */
    SUNMatrix_ID matrix_id() const {return SUNMatGetID(matrix_);};

  private:
    void update_ptrs();

    /**
     * @brief CSC matrix to which all methods are applied
     */
    SUNMatrix matrix_ {nullptr};
    realtype *data_ptr_ {nullptr};
    sunindextype *indexptrs_ptr_ {nullptr};
    sunindextype *indexvals_ptr_ {nullptr};
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
