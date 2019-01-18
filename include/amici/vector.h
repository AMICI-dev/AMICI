#ifndef AMICI_VECTOR_H
#define AMICI_VECTOR_H

#include <vector>

#include <nvector/nvector_serial.h>
#include <amici/exception.h>

namespace amici {

/** AmiVector class provides a generic interface to the NVector_Serial struct */
class AmiVector {
public:
    /**
      * Creates an std::vector<realtype> and attaches the
      * data pointer to a newly created N_Vector_Serial.
      * Using N_VMake_Serial ensures that the N_Vector
      * module does not try to deallocate the data vector
      * when calling N_VDestroy_Serial
      * @param length number of elements in vector
      * @return new AmiVector instance
      */

    explicit AmiVector(const long int length)
        : vec(static_cast<decltype(vec)::size_type>(length), 0.0),
          nvec(N_VMake_Serial(length,vec.data()))
    {
    }

    /** constructor from std::vector, copies data from std::vector
      * and constructs an nvec that points to the data
      * @param rvec vector from which the data will be copied
      * @return new AmiVector instance
      */
    explicit AmiVector(std::vector<realtype> rvec)
        : vec(std::move(rvec)),
          nvec(N_VMake_Serial(static_cast<long int>(vec.size()), vec.data()))
    {
    }

    /** copy constructor
      * @param vold vector from which the data will be copied
      */
    AmiVector(const AmiVector& vold): vec(vold.vec) {
        nvec = N_VMake_Serial(static_cast<long int>(vold.vec.size()), vec.data());
    }

    /** copy-move assignment operator
      * @param other right hand side
      * @return left hand side
      */
    AmiVector& operator=(AmiVector const& other) {
        vec = other.vec;
        if(nvec)
            N_VDestroy_Serial(nvec);
        nvec = N_VMake_Serial(static_cast<long int>(vec.size()),vec.data());
        return *this;
    }

    /** data accessor
      * @return pointer to data array
      */
    realtype *data() {
        return vec.data();
    }

    /** const data accessor
      * @return const pointer to data array
      */
    const realtype *data() const {
        return vec.data();
    }

    /** N_Vector accessor
      * @return N_Vector
      */
    N_Vector getNVector() const {
        return nvec;
    }

    /** Vector accessor
      * @return Vector
      */
    std::vector<realtype> const& getVector() const {
        return vec;
    }

    /** returns the length of the vector
      * @return length
      */
    int getLength() const {
        return static_cast<int>(vec.size());
    }

    /** resets the Vector by filling with zero values
      */
    void reset() {
        set(0.0);
    }

    /** changes the sign of data elements
      */
    void minus() {
        for(std::vector<realtype>::iterator it = vec.begin();
            it != vec.end(); ++it)
         *it = -*it;
    }

    /** sets all data elements to a specific value
      * @param val value for data elements
      */
    void set(realtype val) {
        std::fill(vec.begin(), vec.end(), val);
    }

    /** accessor to data elements of the vector
      * @param pos index of element
      * @return element
      */
    realtype& operator[](int pos) {
        return vec.at(static_cast<decltype(vec)::size_type>(pos));
    }

    /** accessor to data elements of the vector
      * @param pos index of element
      * @return element
      */
    realtype& at(int pos) {
        return vec.at(static_cast<decltype(vec)::size_type>(pos));
    }

    /** accessor to data elements of the vector
     * @param pos index of element
     * @return element
     */
    const realtype& at(int pos) const {
        return vec.at(static_cast<decltype(vec)::size_type>(pos));
    }

    ~AmiVector(){
        N_VDestroy_Serial(nvec);
    }

private:
    /** main data storage */
    std::vector<realtype> vec;
    /** N_Vector, will be synchronised such that it points to
      * data in vec */
    N_Vector nvec = nullptr;
};


/** AmiVectorArray class.
     provides a generic interface to arrays of NVector_Serial structs
*/
class AmiVectorArray {
public:
    /** creates an std::vector<realype> and attaches the
      * data pointer to a newly created N_VectorArray
      * using CloneVectorArrayEmpty ensures that the N_Vector
      * module does not try to deallocate the data vector
      * when calling N_VDestroyVectorArray_Serial
      * @param length_inner length of vectors
      * @param length_outer number of vectors
      * @return New AmiVectorArray instance
      */
    AmiVectorArray(long int length_inner, long int length_outer)
        : vec_array(static_cast<decltype(vec_array)::size_type>(length_outer),
                    AmiVector(length_inner))
    {
        nvec_array = new N_Vector[length_outer];
        for (int idx = 0; idx < length_outer; idx++) {
            nvec_array[idx] = vec_array.at(static_cast<decltype(vec_array)::size_type>(idx)).getNVector();
        }
    }

    /** copy constructor
      * @param vaold object to copy from
      * @return new AmiVectorArray instance
      */
    AmiVectorArray(const AmiVectorArray& vaold) : vec_array(vaold.vec_array) {
        nvec_array = new N_Vector[vaold.getLength()];
        for (int idx = 0; idx < vaold.getLength(); idx++) {
            nvec_array[idx] = vec_array.at(static_cast<decltype(vec_array)::size_type>(idx)).getNVector();
        }
    }

    /** accessor to data of AmiVector elements
      * @param pos index of AmiVector
      * @return pointer to data array
      */
    realtype *data(int pos) {
        return vec_array.at(static_cast<decltype(vec_array)::size_type>(pos)).data();
    }

    /** const accessor to data of AmiVector elements
      * @param pos index of AmiVector
      * @return const pointer to data array
      */
    const realtype *data(int pos) const {
        return vec_array.at(static_cast<decltype(vec_array)::size_type>(pos)).data();
    }

    /** accessor to elements of AmiVector elements
     * @param ipos inner index in AmiVector
     * @param jpos outer index in AmiVectorArray
     * @return element
     */
    realtype& at(int ipos, int jpos) {
        return vec_array.at(static_cast<decltype(vec_array)::size_type>(jpos)).at(ipos);
    }

    /** accessor to elements of AmiVector elements
      * @param ipos inner index in AmiVector
      * @param jpos outer index in AmiVectorArray
      * @return element
      */
    const realtype& at(int ipos, int jpos) const {
        return vec_array.at(static_cast<decltype(vec_array)::size_type>(jpos)).at(ipos);
    }

    /** accessor to NVectorArray
      * @return N_VectorArray
      */
    N_Vector *getNVectorArray() {
        return nvec_array;
    }

    /** accessor to NVector element
      * @param pos index of corresponding AmiVector
      * @return N_Vector
      */
    N_Vector getNVector(int pos) {
        return nvec_array[pos];
    }

    /** accessor to AmiVector elements
      * @param pos index of AmiVector
      * @return AmiVector
      */
    AmiVector& operator[](int pos) {
        return vec_array.at(static_cast<decltype(vec_array)::size_type>(pos));
    }

    /** const accessor to AmiVector elements
      * @param pos index of AmiVector
      * @return const AmiVector
      */
    const AmiVector& operator[](int pos) const {
        return vec_array.at(static_cast<decltype(vec_array)::size_type>(pos));
    }

    /** length of AmiVectorArray
      * @return length
      */
    int getLength() const {
        return static_cast<int>(vec_array.size());
    }

    /** resets every AmiVector in AmiVectorArray */
    void reset() {
        for(auto &v: vec_array)
            v.reset();
    }

    /** flattens the AmiVectorArray to a vector in row-major format
     * @param vec vector into which the AmiVectorArray will be flattened. Must
     * have length equal to number of elements.
     */
    void flatten_to_vector(std::vector<realtype>& vec) const {
        int n_outer = vec_array.size();
        if(n_outer == 0)
            return; //nothing to do ...
        int n_inner = vec_array.at(0).getLength();

        if (static_cast<int>(vec.size()) != n_inner * n_outer) {
            throw AmiException("Dimension of AmiVectorArray (%ix%i) does not "
                               "match target vector dimension (%u)",
                               n_inner, n_outer, vec.size());
        }

        for (int outer = 0; outer < n_outer; ++outer) {
            for (int inner = 0; inner < n_inner; ++inner)
                vec.at(inner + outer * n_inner) = this->at(inner,outer);
        }
    }

    ~AmiVectorArray(){
        delete[] nvec_array;
    }

private:
    /** main data storage */
    std::vector<AmiVector> vec_array;
    /** N_Vector array, will be synchronised such that it points to
      * respective elements in the vec_array
      */
    N_Vector *nvec_array = nullptr;
};

}


#endif /* AMICI_VECTOR_H */
