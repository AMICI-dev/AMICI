//
//  amici_vector.h

#include <vector>
#include <nvector/nvector_serial.h>

#ifndef amici_vector_h
#define amici_vector_h

namespace amici {
    
    /** AmiVector class.
     provides a generic interface to the NVector_Serial struct
     */
    class AmiVector {
    public:
        /** default constructor
         * creates an std::vector<realype> and attaches the 
         * data pointer to a newly created N_Vector_Serial
         * using N_VMake_Serial ensures that the N_Vector
         * module does not try to deallocate the data vector
         * when calling N_VDestroy_Serial
         * @param length number of elements in vector
         * @return new AmiVector instance
         */
        AmiVector(const long int length) : vec(length,0.0) {
            nvec = N_VMake_Serial(length,vec.data());
        };
        
        /** constructor from vector, copies data from vector
         * and constructs an nvec that points to the data
         * @param rvec vector from which the data will be copied
         * @return new AmiVector instance
         */
        AmiVector(std::vector<realtype> rvec) {
            vec = rvec;
            nvec = N_VMake_Serial(rvec.size(),rvec.data());
        };
        
        /** swap functions for copy-and-swap idiom, swaps data between two AmiVectors
         * see https://stackoverflow.com/questions/3279543/what-is-the-copy-and-swap-idiom
         * @param first will be switched with second
         * @param second will be switched with first
         */
        friend void swap( AmiVector& first, AmiVector& second ) {
            std::swap(first.vec,second.vec);
            if(first.nvec)
                N_VDestroy_Serial(first.nvec);
            first.nvec = N_VMake_Serial(first.vec.size(),first.vec.data());
            if(second.nvec)
                N_VDestroy_Serial(second.nvec);
            second.nvec = N_VMake_Serial(second.vec.size(),second.vec.data());
        };
        
        /** copy constructor
         * @param vold vector from which the data will be copied
         */
        AmiVector(const AmiVector& vold) {
            vec = vold.vec;
            nvec = N_VMake_Serial(vold.vec.size(),vec.data());
        };
        
        /** copy-move assignment operator
         * @param other right hand side
         * @return left hand side
         */
        AmiVector& operator=(AmiVector& other) {
            vec = other.vec;
            if(nvec)
                N_VDestroy_Serial(nvec);
            nvec = N_VMake_Serial(vec.size(),vec.data());
            return *this;
        };
        
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
        N_Vector getNVector() {
            return nvec;
        }
        
        /** const N_Vector accessor
         * @return const N_Vector
         */
        const N_Vector getNVector() const {
            return nvec;
        }
        
        /** returns the length of the vector
         * @return length
         */
        const int getLength() const {
            return vec.size();
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
            return vec.at(pos);
        }
        
        /** accessor to data elements of the vector
         * @param pos index of element
         * @return element
         */
        realtype& at(int pos) {
            return vec.at(pos);
        }
        
        /** default destructor
         * creates an std::vector<realype> and attaches the
         * data pointer to a newly created N_Vector_Serial
         */
        ~AmiVector(){
            N_VDestroy_Serial(nvec);
        };
    private:
        /** main data storage
         */
        std::vector<realtype> vec;
        /** N_Vector, will be synchronised such that it points to
         * data in vec
         */
        N_Vector nvec = nullptr;
    };
    
    /** AmiVectorArray class.
     provides a generic interface to arrays of NVector_Serial structs
     */
    class AmiVectorArray {
    public:
        /** default constructor
         * creates an std::vector<realype> and attaches the
         * data pointer to a newly created N_VectorArray
         * using CloneVectorArrayEmpty ensures that the N_Vector
         * module does not try to deallocate the data vector
         * when calling N_VDestroyVectorArray_Serial
         * @param length_inner length of vectors
         * @param length_outer number of vectors
         * @return New AmiVectorArray instance
         */
        AmiVectorArray(const long int length_inner, const long int length_outer)
        : vec_array(length_outer,AmiVector(length_inner))
        {
            nvec_array = new N_Vector[length_outer];
            for (int idx = 0; idx < length_outer; idx++) {
                nvec_array[idx] = vec_array.at(idx).getNVector();
            }
        };
        
        /** copy constructor
         * @param vaold object to copy from
         * @return new AmiVectorArray instance
         */
        AmiVectorArray(const AmiVectorArray& vaold) : vec_array(vaold.vec_array) {
            nvec_array = new N_Vector[vaold.getLength()];
            for (int idx = 0; idx < vaold.getLength(); idx++) {
                nvec_array[idx] = vec_array.at(idx).getNVector();
            }
        }
        
        /** accessor to data of AmiVector elements
         * @param pos index of AmiVector
         * @return pointer to data array
         */
        realtype *data(int pos) {
            return vec_array.at(pos).data();
        }
        
        /** const accessor to data of AmiVector elements
         * @param pos index of AmiVector
         * @return const pointer to data array
         */
        const realtype *data(int pos) const {
            return vec_array.at(pos).data();
        }
        
        /** accessor to elements of AmiVector elements
         * @param ipos inner index in AmiVector
         * @param jpos outer index in AmiVectorArray
         * @return element
         */
        realtype& at(int ipos, int jpos) {
            return vec_array.at(jpos)[ipos];
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
            return vec_array.at(pos);
        }
        
        /** const accessor to AmiVector elements
         * @param pos index of AmiVector
         * @return const AmiVector
         */
        const AmiVector& operator[](int pos) const {
            return vec_array.at(pos);
        }
        
        /** length of AmiVectorArray
         * @return length
         */
        const int getLength() const {
            return vec_array.size();
        }
        
        /** resets every AmiVector in AmiVectorArray
         */
        void reset() {
            for(std::vector<AmiVector>::iterator it = vec_array.begin();
                it != vec_array.end(); ++it)
                it->reset();
        }
        
        /** default destructor
         */
        ~AmiVectorArray(){
            delete[] nvec_array;
        }
    private:
        /** main data storage
         */
        std::vector<AmiVector> vec_array;
        /** N_Vector array, will be synchronised such that it points to
         * respective elements in the vec_array
         */
        N_Vector *nvec_array = nullptr;
    };
    
}


#endif /* amici_nvector_h */
