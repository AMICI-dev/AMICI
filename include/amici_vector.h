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
        
        realtype *data() {
            return vec.data();
        }
        
        const realtype *data() const {
            return vec.data();
        }
        
        N_Vector getNVector() {
            return nvec;
        }
        
        int getLength() {
            return vec.size();
        }
        
        void reset() {
            std::fill(vec.begin(), vec.end(), 0.0);
        }
        
        realtype& operator[](int pos) {
            return vec.at(pos);
        }
        
        /*! default destructor
         * creates an std::vector<realype> and attaches the
         * data pointer to a newly created N_Vector_Serial
         */
        virtual ~AmiVector(){
            N_VDestroy_Serial(nvec);
        };
    private:
        std::vector<realtype> vec;
        N_Vector nvec;
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
        : vec_array(length_outer,std::vector<realtype>(length_inner,0.0))
        {
            N_Vector nvec = N_VMake_Serial(length_inner,vec_array.at(0).data());
            nvec_array = N_VCloneVectorArrayEmpty_Serial(length_outer,nvec);
            for (int idx = 0; idx < length_outer; idx++)
                NV_DATA_S(nvec_array[idx]) = vec_array.at(idx).data();
        };
        
        realtype *data(int idx) {
            return vec_array.at(idx).data();
        }
        
        realtype& at(int ipos, int jpos) {
            return vec_array.at(jpos).at(ipos);
        }
        
        const realtype *data(int idx) const {
            return vec_array.at(idx).data();
        }
        
        N_Vector *getNVectorArray() {
            return nvec_array;
        }
        
        int getLength() {
            return vec_array.size();
        }
        
        void reset() {
            for(std::vector<std::vector<realtype>>::iterator it = vec_array.begin();
                it != vec_array.end(); ++it)
                std::fill(it->begin(), it->end(), 0.0);
        }
        
        /*! default destructor
         * creates an std::vector<realype> and attaches the
         * data pointer to a newly created N_Vector_Serial
         */
        virtual ~AmiVectorArray(){
            N_VDestroyVectorArray_Serial(nvec_array,vec_array.size());
        };
    private:
        std::vector<std::vector<realtype>> vec_array;
        N_Vector *nvec_array;
    };
    
}


#endif /* amici_nvector_h */
