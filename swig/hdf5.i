%module hdf5

// Add necessary symbols to generated header
%{
#include "amici_hdf5.h"
using namespace amici;
%}

%include std_unique_ptr.i

wrap_unique_ptr(ExpDataPtr, amici::ExpData)

%rename(readModelDataFromHDF5) readModelDataFromHDF5Swig;
%inline %{

    Model_model_events *getModelSwig() {
        return new Model_model_events();
    }
    // Helper function to dereference pointers within python
    template <typename T>
    T& dereference(T* ptr)
    {
        return *ptr;
    }
%}
%ignore getModel;


// Process symbols in header

%include "amici_hdf5.h"
