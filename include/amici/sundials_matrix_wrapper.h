#ifndef AMICI_SUNDIALS_MATRIX_WRAPPER_H
#define AMICI_SUNDIALS_MATRIX_WRAPPER_H

#include <sunmatrix/sunmatrix_band.h>   // SUNMatrix_Band
#include <sunmatrix/sunmatrix_dense.h>  // SUNMatrix_Dense
#include <sunmatrix/sunmatrix_sparse.h> // SUNMatrix_Sparse

#include <nvector/nvector_serial.h>

#include <vector>

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
    SUNMatrixWrapper(SUNMatrixWrapper &&other) noexcept;

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
    SUNMatrixWrapper &operator=(SUNMatrixWrapper &&other) noexcept;

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
     * @brief std::vector interface for multiply
     * @param c output vector, may already contain values
     * @param b multiplication vector
     */
    void multiply(std::vector<realtype> &c, const std::vector<realtype> &b) const;

    /**
     * @brief N_Vector interface for multiply
     * @param c output vector, may already contain values
     * @param b multiplication vector
     */
    void multiply(N_Vector c, const N_Vector b) const;

    /**
     * @brief Perform matrix vector multiplication c += A*b
     * @param c output vector, may already contain values
     * @param b multiplication vector
     */
    void multiply(realtype *c, const realtype *b) const;

  private:
    void update_ptrs();
    
    SUNMatrix matrix = nullptr;
    realtype *data_ptr = nullptr;
    sunindextype *indexptrs_ptr = nullptr;
    sunindextype *indexvals_ptr = nullptr;
};

} // namespace amici

#endif // AMICI_SUNDIALS_MATRIX_WRAPPER_H
