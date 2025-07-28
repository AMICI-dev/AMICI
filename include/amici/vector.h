#ifndef AMICI_VECTOR_H
#define AMICI_VECTOR_H

#include <ostream>
#include <type_traits>
#include <vector>

#include <amici/exception.h>

#include <nvector/nvector_serial.h>

#include <gsl/gsl-lite.hpp>
#include <sundials/sundials_context.hpp>
namespace amici {
class AmiVector;
}

// for serialization friend
namespace boost::serialization {
template <class Archive>
void serialize(Archive& ar, amici::AmiVector& v, unsigned int version);
} // namespace boost::serialization

namespace amici {

/** Since const N_Vector is not what we want */
using const_N_Vector = std::add_const_t<std::remove_pointer_t<N_Vector>>*;

inline realtype const* N_VGetArrayPointerConst(const_N_Vector x) {
    return N_VGetArrayPointer(const_cast<N_Vector>(x));
}

/** AmiVector class provides a generic interface to the NVector_Serial struct */
class AmiVector {
  public:
    /**
     * @brief Default constructor
     */
    AmiVector() = default;

    /**
     * @brief Construct zero-initialized vector of the given size.
     *
     * Creates an std::vector<realtype> and attaches the
     * data pointer to a newly created N_Vector_Serial.
     * Using N_VMake_Serial ensures that the N_Vector
     * module does not try to deallocate the data vector
     * when calling N_VDestroy_Serial
     * @param length number of elements in vector
     * @param sunctx SUNDIALS context
     */
    explicit AmiVector(long int const length, SUNContext const sunctx)
        : vec_(static_cast<decltype(vec_)::size_type>(length), 0.0)
        , nvec_(N_VMake_Serial(length, vec_.data(), sunctx)) {}

    /**
     * @brief Constructor from std::vector
     *
     * Moves data from std::vector and constructs an nvec that points to the
     * data
     * @param rvec vector from which the data will be moved
     * @param sunctx SUNDIALS context
     */
    explicit AmiVector(std::vector<realtype> rvec, SUNContext const sunctx)
        : vec_(std::move(rvec))
        , nvec_(N_VMake_Serial(
              gsl::narrow<long int>(vec_.size()), vec_.data(), sunctx
          )) {}

    /** Copy data from gsl::span and constructs a vector
     * @brief constructor from gsl::span,
     * @param rvec vector from which the data will be copied
     * @param sunctx SUNDIALS context
     */
    explicit AmiVector(
        gsl::span<realtype const> const rvec, SUNContext const sunctx
    )
        : AmiVector(std::vector(rvec.begin(), rvec.end()), sunctx) {}

    /**
     * @brief copy constructor
     * @param vold vector from which the data will be copied
     */
    AmiVector(AmiVector const& vold)
        : vec_(vold.vec_) {
        if (vold.nvec_ == nullptr) {
            nvec_ = nullptr;
            return;
        }
        nvec_ = N_VMake_Serial(
            gsl::narrow<long int>(vec_.size()), vec_.data(), vold.nvec_->sunctx
        );
    }

    /**
     * @brief Copy from N_Vector
     * @param vold vector from which the data will be copied
     */
    AmiVector(N_Vector const& vold)
        : AmiVector(
              gsl::span(
                  N_VGetArrayPointerConst(vold), N_VGetLength_Serial(vold)
              ),
              vold->sunctx
          ) {}

    /**
     * @brief move constructor
     * @param other vector from which the data will be moved
     */
    AmiVector(AmiVector&& other) noexcept
        : vec_(std::move(other.vec_)) {
        synchroniseNVector(other.get_ctx());
    }

    /**
     * @brief destructor
     */
    ~AmiVector();

    /**
     * @brief copy assignment operator
     * @param other right hand side
     * @return left hand side
     */
    AmiVector& operator=(AmiVector const& other);

    /**
     * @brief operator *= (element-wise multiplication)
     * @param multiplier multiplier
     * @return result
     */
    AmiVector& operator*=(AmiVector const& multiplier) {
        N_VProd(
            getNVector(), const_cast<N_Vector>(multiplier.getNVector()),
            getNVector()
        );
        return *this;
    }

    /**
     * @brief operator /= (element-wise division)
     * @param divisor divisor
     * @return result
     */
    AmiVector& operator/=(AmiVector const& divisor) {
        N_VDiv(
            getNVector(), const_cast<N_Vector>(divisor.getNVector()),
            getNVector()
        );
        return *this;
    }

    /**
     * @brief Returns an iterator that points to the first element of the
     * vector.
     * @return iterator that points to the first element
     */
    auto begin() { return vec_.begin(); }

    /**
     * @brief Returns an iterator that points to one element after the last
     * element of the vector.
     * @return iterator that points to one element after the last element
     */
    auto end() { return vec_.end(); }

    /**
     * @brief data accessor
     * @return pointer to data array
     */
    realtype* data();

    /**
     * @brief const data accessor
     * @return const pointer to data array
     */
    [[nodiscard]] realtype const* data() const;

    /**
     * @brief N_Vector accessor
     * @return N_Vector
     */
    N_Vector getNVector();

    /**
     * @brief N_Vector accessor
     * @return N_Vector
     */
    [[nodiscard]] const_N_Vector getNVector() const;

    /**
     * @brief Vector accessor
     * @return Vector
     */
    [[nodiscard]] std::vector<realtype> const& getVector() const;

    /**
     * @brief returns the length of the vector
     * @return length
     */
    [[nodiscard]] int getLength() const;

    /**
     * @brief fills vector with zero values
     */
    void zero();

    /**
     * @brief changes the sign of data elements
     */
    void minus();

    /**
     * @brief sets all data elements to a specific value
     * @param val value for data elements
     */
    void set(realtype val);

    /**
     * @brief accessor to data elements of the vector
     * @param pos index of element
     * @return element
     */
    realtype& operator[](int pos);

    /**
     * @brief accessor to data elements of the vector
     * @param pos index of element
     * @return element
     */
    realtype const& operator[](int pos) const;

    /**
     * @brief accessor to data elements of the vector
     * @param pos index of element
     * @return element
     */
    realtype& at(int pos);

    /**
     * @brief accessor to data elements of the vector
     * @param pos index of element
     * @return element
     */
    [[nodiscard]] realtype const& at(int pos) const;

    /**
     * @brief copies data from another AmiVector
     * @param other data source
     */
    void copy(AmiVector const& other);

    /**
     * @brief Take absolute value (in-place)
     */
    void abs() { N_VAbs(getNVector(), getNVector()); }

    /**
     * @brief Serialize AmiVector (see boost::serialization::serialize)
     * @param ar Archive to serialize to
     * @param v Data to serialize
     * @param version Version number
     */
    template <class Archive>
    friend void boost::serialization::serialize(
        Archive& ar, AmiVector& v, unsigned int version
    );

    /**
     * @brief Get SUNContext
     * @return The current SUNContext or nullptr, if this AmiVector is empty
     */
    [[nodiscard]] SUNContext get_ctx() const {
        return nvec_ == nullptr ? nullptr : nvec_->sunctx;
    }

    /**
     * @brief Set SUNContext
     *
     * If this AmiVector is non-empty, changes the current SUNContext of the
     * associated N_Vector. If empty, do nothing.
     *
     * @param ctx SUNDIALS context to set
     */
    void set_ctx(SUNContext ctx) {
        if (nvec_)
            nvec_->sunctx = ctx;
    }

  private:
    /** main data storage */
    std::vector<realtype> vec_;

    /** N_Vector, will be synchronized such that it points to data in vec */
    N_Vector nvec_{nullptr};

    /**
     * @brief reconstructs nvec such that data pointer points to vec data array
     * @param sunctx SUNDIALS context
     */
    void synchroniseNVector(SUNContext sunctx);
};

/**
 * @brief Output formatter for AmiVector.
 * @param os output stream
 * @param v AmiVector to output
 * @return os
 */
inline std::ostream& operator<<(std::ostream& os, AmiVector const& v) {
    os << "[";
    for (int i = 0; i < v.getLength(); ++i) {
        if (i > 0)
            os << ", ";
        os << v.at(i);
    }
    os << "]";
    return os;
}

/**
 * @brief AmiVectorArray class.
 *
 * Provides a generic interface to arrays of NVector_Serial structs
 */
class AmiVectorArray {
  public:
    /**
     * @brief Default constructor
     */
    AmiVectorArray() = default;

    /**
     * Creates an std::vector<realype> and attaches the
     * data pointer to a newly created N_VectorArray
     * using CloneVectorArrayEmpty ensures that the N_Vector
     * module does not try to deallocate the data vector
     * when calling N_VDestroyVectorArray_Serial
     * @brief empty constructor
     * @param length_inner length of vectors
     * @param length_outer number of vectors
     * @param sunctx SUNDIALS context
     */
    AmiVectorArray(
        long int length_inner, long int length_outer, SUNContext sunctx
    );

    /**
     * @brief copy constructor
     * @param vaold object to copy from
     */
    AmiVectorArray(AmiVectorArray const& vaold);

    ~AmiVectorArray() = default;

    /**
     * @brief copy assignment operator
     * @param other right hand side
     * @return left hand side
     */
    AmiVectorArray& operator=(AmiVectorArray const& other);

    /**
     * @brief accessor to data of AmiVector elements
     * @param pos index of AmiVector
     * @return pointer to data array
     */
    realtype* data(int pos);

    /**
     * @brief const accessor to data of AmiVector elements
     * @param pos index of AmiVector
     * @return const pointer to data array
     */
    [[nodiscard]] realtype const* data(int pos) const;

    /**
     * @brief accessor to elements of AmiVector elements
     * @param ipos inner index in AmiVector
     * @param jpos outer index in AmiVectorArray
     * @return element
     */
    realtype& at(int ipos, int jpos);

    /**
     * @brief const accessor to elements of AmiVector elements
     * @param ipos inner index in AmiVector
     * @param jpos outer index in AmiVectorArray
     * @return element
     */
    [[nodiscard]] realtype const& at(int ipos, int jpos) const;

    /**
     * @brief accessor to NVectorArray
     * @return N_VectorArray
     */
    N_Vector* getNVectorArray();

    /**
     * @brief accessor to NVector element
     * @param pos index of corresponding AmiVector
     * @return N_Vector
     */
    N_Vector getNVector(int pos);

    /**
     * @brief const accessor to NVector element
     * @param pos index of corresponding AmiVector
     * @return N_Vector
     */
    [[nodiscard]] const_N_Vector getNVector(int pos) const;

    /**
     * @brief accessor to AmiVector elements
     * @param pos index of AmiVector
     * @return AmiVector
     */
    AmiVector& operator[](int pos);

    /**
     * @brief const accessor to AmiVector elements
     * @param pos index of AmiVector
     * @return const AmiVector
     */
    AmiVector const& operator[](int pos) const;

    /**
     * @brief length of AmiVectorArray
     * @return length
     */
    int getLength() const;

    /**
     * @brief set every AmiVector in AmiVectorArray to zero
     */
    void zero();

    /**
     * @brief flattens the AmiVectorArray to a vector in row-major format
     * @param vec vector into which the AmiVectorArray will be flattened. Must
     * have length equal to number of elements.
     */
    void flatten_to_vector(std::vector<realtype>& vec) const;

    /**
     * @brief copies data from another AmiVectorArray
     * @param other data source
     */
    void copy(AmiVectorArray const& other);

    /**
     * @brief Set SUNContext
     *
     * If any AmiVector is non-empty, this changes the current SUNContext of the
     * associated N_Vector. If empty, do nothing.
     *
     * @param ctx SUNDIALS context to set
     */
    void set_ctx(SUNContext ctx) {
        for (auto& vec : vec_array_) {
            vec.set_ctx(ctx);
        }
    }

  private:
    /** main data storage */
    std::vector<AmiVector> vec_array_;

    /**
     * N_Vector array, will be synchronized such that it points to
     * respective elements in the vec_array
     */
    std::vector<N_Vector> nvec_array_;
};

/**
 * @brief Output formatter for AmiVectorArray.
 * @param os output stream
 * @param arr AmiVectorArray to output
 * @return os
 */
inline std::ostream& operator<<(std::ostream& os, AmiVectorArray const& arr) {
    os << "[";
    for (int i = 0; i < arr.getLength(); ++i) {
        if (i > 0)
            os << ", ";
        os << arr[i];
    }
    os << "]";
    return os;
}

/**
 * @brief Computes z = a*x + b*y
 * @param a coefficient for x
 * @param x a vector
 * @param b coefficient for y
 * @param y another vector with same size as x
 * @param z result vector of same size as x and y
 */
inline void linearSum(
    realtype const a, AmiVector const& x, realtype const b, AmiVector const& y,
    AmiVector& z
) {
    N_VLinearSum(
        a, const_cast<N_Vector>(x.getNVector()), b,
        const_cast<N_Vector>(y.getNVector()), z.getNVector()
    );
}

/**
 * @brief Compute dot product of x and y
 * @param x vector
 * @param y vector
 * @return dot product of x and y
 */
inline realtype dotProd(AmiVector const& x, AmiVector const& y) {
    return N_VDotProd(
        const_cast<N_Vector>(x.getNVector()),
        const_cast<N_Vector>(y.getNVector())
    );
}

} // namespace amici

namespace gsl {
/**
 * @brief Create span from N_Vector
 * @param nv
 * @return
 */
inline span<amici::realtype> make_span(N_Vector nv) {
    return span(N_VGetArrayPointer(nv), N_VGetLength_Serial(nv));
}

/**
 * @brief Create span from AmiVector
 * @param av
 *
 */
inline span<amici::realtype const> make_span(amici::AmiVector const& av) {
    return make_span(av.getVector());
}

} // namespace gsl

#endif /* AMICI_VECTOR_H */
