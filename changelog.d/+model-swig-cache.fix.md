Rebuilding a model extension in an existing model directory with a
different SWIG version now reruns SWIG. Previously, the CMake cache
could keep the library directory and version of the previously used
SWIG, and the SWIG-generated code was not updated.
