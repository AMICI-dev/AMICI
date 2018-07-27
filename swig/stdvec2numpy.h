namespace amici {

/**
 * @brief Convert 1D array to *non-owning* numpy ndarray.
 * @param vec
 * @param dim1
 * @return
 */
PyObject* stdVec2ndarray(std::vector<double>& vec, int dim1) {
    if (vec.size() != (unsigned) dim1) throw std::runtime_error("Size mismatch in stdVec2ndarray");
    npy_intp dims[1] = { dim1 };
    PyObject * array = PyArray_SimpleNewFromData(1, dims,  NPY_DOUBLE, vec.data());
    if (!array) throw std::runtime_error("Unknown failure in stdVec2ndarray");;
    return array;
}

/**
 * @brief Convert row-major flattened 2D array to *non-owning* numpy ndarray.
 * @param vec
 * @param dim1
 * @param dim2
 * @return
 */
PyObject* stdVec2ndarray(std::vector<double>& vec, int dim1, int dim2) {
    if (vec.size() != (unsigned) dim1 * dim2) throw std::runtime_error("Size mismatch in stdVec2ndarray");
    npy_intp dims[2] = { dim1, dim2 };
    PyObject * array = PyArray_SimpleNewFromData(2, dims,  NPY_DOUBLE, vec.data());
    if (!array) throw std::runtime_error("Unknown failure in stdVec2ndarray");;
    return array;
}

/**
 * @brief Convert row-major flattened 3D array to *non-owning* numpy ndarray.
 * @param vec
 * @param dim1
 * @param dim2
 * @param dim3
 * @return
 */
PyObject* stdVec2ndarray(std::vector<double>& vec, int dim1, int dim2, int dim3) {
    if (vec.size() != (unsigned) dim1 * dim2 * dim3) throw std::runtime_error("Size mismatch in stdVec2ndarray");
    npy_intp dims[3] = { dim1, dim2, dim3 };
    PyObject * array = PyArray_SimpleNewFromData(3, dims,  NPY_DOUBLE, vec.data());
    if (!array) throw std::runtime_error("Unknown failure in stdVec2ndarray");;
    return array;
}

/**
 * @brief Convert row-major flattened 2D array to *non-owning* numpy ndarray.
 * @param vec
 * @param dim1
 * @param dim2
 * @param dim3
 * @param dim4
 * @return
 */
PyObject* stdVec2ndarray(std::vector<double>& vec, int dim1, int dim2, int dim3, int dim4) {
    if (vec.size() != (unsigned) dim1 * dim2 * dim3 * dim4) throw std::runtime_error("Size mismatch in stdVec2ndarray");
    npy_intp dims[4] = { dim1, dim2, dim3, dim4 };
    PyObject * array = PyArray_SimpleNewFromData(4, dims,  NPY_DOUBLE, vec.data());
    if (!array) throw std::runtime_error("Unknown failure in stdVec2ndarray");;
    return array;
}


/**
 * @brief Convert 1D array to *non-owning* numpy ndarray.
 * @param vec
 * @param dim1
 * @return
 */
PyObject* stdVec2ndarray(std::vector<int>& vec, int dim1) {
    if (vec.size() != (unsigned) dim1) throw std::runtime_error("Size mismatch in stdVec2ndarray");
    npy_intp dims[1] = { dim1 };
    PyObject * array = PyArray_SimpleNewFromData(1, dims,  NPY_INT, vec.data());
    if (!array) throw std::runtime_error("Unknown failure in stdVec2ndarray");;
    return array;
}

/**
 * @brief Convert row-major flattened 2D array to *non-owning* numpy ndarray.
 * @param vec
 * @param dim1
 * @param dim2
 * @return
 */
PyObject* stdVec2ndarray(std::vector<int>& vec, int dim1, int dim2) {
    if (vec.size() != (unsigned) dim1 * dim2) throw std::runtime_error("Size mismatch in stdVec2ndarray");
    npy_intp dims[2] = { dim1, dim2 };
    PyObject * array = PyArray_SimpleNewFromData(2, dims,  NPY_INT, vec.data());
    if (!array) throw std::runtime_error("Unknown failure in stdVec2ndarray");;
    return array;
}

/**
 * @brief Convert row-major flattened 3D array to *non-owning* numpy ndarray.
 * @param vec
 * @param dim1
 * @param dim2
 * @param dim3
 * @return
 */
PyObject* stdVec2ndarray(std::vector<int>& vec, int dim1, int dim2, int dim3) {
    if (vec.size() != (unsigned) dim1 * dim2 * dim3) throw std::runtime_error("Size mismatch in stdVec2ndarray");
    npy_intp dims[3] = { dim1, dim2, dim3 };
    PyObject * array = PyArray_SimpleNewFromData(3, dims,  NPY_INT, vec.data());
    if (!array) throw std::runtime_error("Unknown failure in stdVec2ndarray");;
    return array;
}

/**
 * @brief Convert row-major flattened 2D array to *non-owning* numpy ndarray.
 * @param vec
 * @param dim1
 * @param dim2
 * @param dim3
 * @param dim4
 * @return
 */
PyObject* stdVec2ndarray(std::vector<int>& vec, int dim1, int dim2, int dim3, int dim4) {
    if (vec.size() != (unsigned) dim1 * dim2 * dim3 * dim4) throw std::runtime_error("Size mismatch in stdVec2ndarray");
    npy_intp dims[4] = { dim1, dim2, dim3, dim4 };
    PyObject * array = PyArray_SimpleNewFromData(4, dims,  NPY_INT, vec.data());
    if (!array) throw std::runtime_error("Unknown failure in stdVec2ndarray");;
    return array;
}

}
