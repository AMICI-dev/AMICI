cd .\SuiteSparse\SuiteSparse_config

mingw64-make library

cd .\..\AMD

mingw64-make library

cd .\..\BTF

mingw64-make library

cd .\..\CAMD

mingw64-make library

cd .\..\COLAMD

mingw64-make library

cd .\..\KLU

mingw64-make library


mkdir .\..\..\sundials\build\
cd .\..\..\sundials\build\

cmake -DCMAKE_INSTALL_PREFIX=".\..\..\build\sundials" \
-DBUILD_ARKODE=OFF \
-DBUILD_CVODE=OFF \
-DBUILD_IDA=OFF \
-DBUILD_KINSOL=OFF \
-DBUILD_SHARED_LIBS=ON \
-DBUILD_STATIC_LIBS=OFF \
-DEXAMPLES_ENABLE=OFF \
-DEXAMPLES_INSTALL=OFF \
-DKLU_ENABLE=ON \
-DKLU_LIBRARY_DIR="$.\..\..\SuiteSparse\lib" \
-DKLU_INCLUDE_DIR="$.\..\..\SuiteSparse\include" \
.. 

make 
make install

cd ..\..\

cmake CMakeLists.txt
make
#cmake -f .\..\models\model_dirac\CMakeLists.txt