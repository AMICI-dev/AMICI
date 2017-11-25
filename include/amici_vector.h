//
//  amici_vector.h

#include <vector>
#include <nvector/nvector_serial.h>

#ifndef amici_vector_h
#define amici_vector_h

namespace amici {
    
    //!  AmiVector class.
    /*!
     provides a generic interface to the NVector_Serial class
     may be extended to other class if necessary
     */
    class AmiVector {
    public:
        /*! default constructor
         * creates an std::vector<realype> and attaches the 
         * data pointer to a newly created N_Vector_Serial
         * using N_VMake_Serial ensures that the N_Vector
         * module does not try to deallocate the data vector
         * when calling N_VDestroy_Serial
         */
        AmiVector(const long int length) : vec(length,0.0) {
            nvec = N_VMake_Serial(length,vec.data());
        };
        
        AmiVector(std::vector<realtype> rvec) {
            vec = rvec;
            nvec = N_VMake_Serial(rvec.size(),rvec.data());
        };
        
        friend void swap( AmiVector& first, AmiVector& second ) {
            std::swap(first.vec,second.vec);
            if(first.nvec)
                N_VDestroy_Serial(first.nvec);
            first.nvec = N_VMake_Serial(first.vec.size(),first.vec.data());
            if(second.nvec)
                N_VDestroy_Serial(second.nvec);
            second.nvec = N_VMake_Serial(second.vec.size(),second.vec.data());
        };
        
        
        AmiVector(const AmiVector& vold) {
            vec = vold.vec;
            nvec = N_VMake_Serial(vold.vec.size(),vec.data());
        };
        
        
        AmiVector& operator=(AmiVector& other) {
            swap(*this,other);
            return *this;
        };
        
        
        realtype *data() {
            return vec.data();
        }
        
        const realtype *data() const {
            return vec.data();
        }
        
        N_Vector getNVector() {
            return nvec;
        }
        
        const N_Vector getNVector() const {
            return nvec;
        }
        
        const int getLength() const {
            return vec.size();
        }
        
        void reset() {
            std::fill(vec.begin(), vec.end(), 0.0);
        }
        
        void minus() {
            for(std::vector<realtype>::iterator it = vec.begin();
                it != vec.end(); ++it)
                *it = -*it;
        }
        
        
        void set(realtype val) {
            std::fill(vec.begin(), vec.end(), val);
        }
        
        realtype& operator[](int pos) {
            return vec.at(pos);
        }
        
        /*! default destructor
         * creates an std::vector<realype> and attaches the
         * data pointer to a newly created N_Vector_Serial
         */
        ~AmiVector(){
            N_VDestroy_Serial(nvec);
        };
    private:
        std::vector<realtype> vec;
        N_Vector nvec = nullptr;
        
        friend class AmiVectorArray;
    };
    
    class AmiVectorArray {
    public:
        /*! default constructor
         * creates an std::vector<realype> and attaches the
         * data pointer to a newly created N_VectorArray
         * using CloneVectorArrayEmpty ensures that the N_Vector
         * module does not try to deallocate the data vector
         * when calling N_VDestroyVectorArray_Serial
         */
        AmiVectorArray(const long int length_inner, const long int length_outer)
        : vec_array(length_outer,AmiVector(length_inner))
        {
            nvec_array = new N_Vector[length_outer];
            for (int idx = 0; idx < length_outer; idx++) {
                nvec_array[idx] = vec_array.at(idx).nvec;
            }
        };
        
        
        AmiVectorArray(const AmiVectorArray& vaold) : vec_array(vaold.vec_array) {
            nvec_array = new N_Vector[vaold.getLength()];
            for (int idx = 0; idx < vaold.getLength(); idx++) {
                nvec_array[idx] = vec_array.at(idx).nvec;
            }
        }
        
        realtype *data(int idx) {
            return vec_array.at(idx).data();
        }
        
        realtype& at(int ipos, int jpos) {
            return vec_array.at(jpos)[ipos];
        }
        
        const realtype *data(int idx) const {
            return vec_array.at(idx).data();
        }
        
        N_Vector *getNVectorArray() {
            return nvec_array;
        }
        
        N_Vector getNVector(int idx) {
            return nvec_array[idx];
        }
        
        AmiVector& operator[](int pos) {
            return vec_array.at(pos);
        }
        
        const AmiVector& operator[](int pos) const {
            return vec_array.at(pos);
        }
        
        
        const int getLength() const {
            return vec_array.size();
        }
        
        void reset() {
            for(std::vector<AmiVector>::iterator it = vec_array.begin();
                it != vec_array.end(); ++it)
                it->reset();
        }
        ~AmiVectorArray(){
            delete[] nvec_array;
        }
    private:
        std::vector<AmiVector> vec_array;
        N_Vector *nvec_array = nullptr;
    };
    
}


#endif /* amici_nvector_h */
