//
//  amici_vector.h


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
         * when calling N_VDestroyVectorArray_Serial
         */
        AmiVector(const long int length) : vec(length,0.0) {
            nvec = N_VMake_Serial(length,vec.data());
        };
        
        realtype *getData() {
            return vec.data;
        }
        
        N_Vector getNvector() {
            return nvec;
        }
        
        int getLength() {
            return vec.size;
        }
        
        void reset() {
            std::fill(vec.begin(), vec.end(), 0.0);
        }
        
        reference operator[]( size_type pos) {
            return vec.at(pos);
        }
        
        /*! default destructor
         * creates an std::vector<realype> and attaches the
         * data pointer to a newly created N_Vector_Serial
         */
        virtual ~AmiVector(){
            N_VDestroyVectorArray_Serial(nvec);
        };
    private:
        std::vector<realtype> vec;
        N_Vector nvec;
    }
    
    class AmiVectorArray {
    public:
        /*! default constructor
         * creates an std::vector<realype> and attaches the
         * data pointer to a newly created N_Vector_Serial
         * using N_VMake_Serial ensures that the N_Vector
         * module does not try to deallocate the data vector
         * when calling N_VDestroyVectorArray_Serial
         */
        AmiVectorArray(const long int length){
            vec.resize(length,0);
            nvec = N_VMake_Serial(length,vec.data());
        };
        
        realtype *getData() {
            return vec.data;
        }
        
        N_Vector getNvectorArray() {
            return nvec_array;
        }
        
        int getLength() {
            return vec.size;
        }
        
        /*! default destructor
         * creates an std::vector<realype> and attaches the
         * data pointer to a newly created N_Vector_Serial
         */
        virtual ~AmiVector(){
            N_VDestroyVectorArray_Serial(nvec);
        };
    private:
        std::vector<std::vector<realtype>> vec_array;
        N_VectorArray nvec_array;
    }
    
}


#endif /* amici_nvector_h */
