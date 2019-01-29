#ifndef AMICI_SUNDIALS_MATRIX_WRAPPER_H
#define AMICI_SUNDIALS_MATRIX_WRAPPER_H

#include <sunmatrix/sunmatrix_band.h>   // SUNMatrix_Band
#include <sunmatrix/sunmatrix_dense.h>  // SUNMatrix_Dense
#include <sunmatrix/sunmatrix_sparse.h> // SUNMatrix_Sparse

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

  private:
    SUNMatrix matrix = nullptr;
};

} // namespace amici

#endif // AMICI_SUNDIALS_MATRIX_WRAPPER_H
