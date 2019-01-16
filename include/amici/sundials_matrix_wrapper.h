#ifndef AMICI_SUNDIALS_MATRIX_WRAPPER_H
#define AMICI_SUNDIALS_MATRIX_WRAPPER_H

#include <sundials/sundials_direct.h> // DlsMat
#include <sundials/sundials_sparse.h> // SlsMat

namespace amici {

/**
 * @brief A RAII wrapper for Sundials SlsMat sparse matrices.
 */
class SlsMatWrapper {
public:
    SlsMatWrapper() = default;

    /**
     * @brief See SparseNewMat in sundials_sparse.h
     * @param M Number of rows
     * @param N Number of columns
     * @param NNZ Number of nonzeros
     * @param sparsetype Sparse type
     */
    SlsMatWrapper(int M, int N, int NNZ, int sparsetype);

    /**
     * @brief SlsMatWrapper
     * @param mat
     */
    explicit SlsMatWrapper(SlsMat mat);

    ~SlsMatWrapper();

    /**
     * @brief Copy constructor
     * @param other
     */
    SlsMatWrapper(const SlsMatWrapper& other);

    /**
     * @brief Move constructor
     * @param other
     */
    SlsMatWrapper(SlsMatWrapper&& other) noexcept;

    /**
     * @brief Copy assignment
     * @param other
     * @return
     */
    SlsMatWrapper& operator=(const SlsMatWrapper& other);

    /**
     * @brief Move assignment
     * @param other
     * @return
     */
    SlsMatWrapper& operator=(SlsMatWrapper&& other) noexcept;

    /**
     * @brief Access raw data
     * @return raw data pointer
     */
    realtype *data();

    /**
     * @brief Get the wrapped SlsMat
     * @return SlsMat
     */
    SlsMat slsmat() const;

private:
    SlsMat matrix = nullptr;
};


/**
 * @brief A RAII wrapper for Sundials DlsMat dense matrices.
 */
class DlsMatWrapper {
public:
    DlsMatWrapper() = default;

    /**
     * @brief See NewDenseMat in sundials_direct.h
     * @param M Number of rows
     * @param N Number of columns
     */
    DlsMatWrapper(long int M, long int N);

    /**
     * @brief DlsMatWrapper
     * @param mat
     */
    explicit DlsMatWrapper(DlsMat mat);

    ~DlsMatWrapper();

    /**
     * @brief Copy constructor
     * @param other
     */
    DlsMatWrapper(const DlsMatWrapper& other);

    /**
     * @brief Move constructor
     * @param other
     */
    DlsMatWrapper(DlsMatWrapper&& other) noexcept;

    /**
     * @brief Copy assignment
     * @param other
     * @return
     */
    DlsMatWrapper& operator=(const DlsMatWrapper& other);

    /**
     * @brief Move assignment
     * @param other
     * @return
     */
    DlsMatWrapper& operator=(DlsMatWrapper&& other) noexcept;

    /**
     * @brief Access raw data
     * @return raw data pointer
     */
    realtype *data();

    /**
     * @brief Get the wrapped DlsMat
     * @return DlsMat
     */
    DlsMat dlsmat() const;

private:
    DlsMat matrix = nullptr;
};

} // namespace amici

#endif // AMICI_SUNDIALS_MATRIX_WRAPPER_H
