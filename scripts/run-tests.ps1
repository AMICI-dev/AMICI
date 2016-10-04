cd .\SuiteSparse\SuiteSparse_config

mingw32-make library -e CC=gcc

cd .\..\AMD

mingw32-make library -e CC=gcc

cd .\..\BTF

mingw32-make library -e CC=gcc

cd .\..\CAMD

mingw32-make library -e CC=gcc

cd .\..\COLAMD

mingw32-make library -e CC=gcc

cd .\..\KLU

mingw32-make library -e CC=gcc


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

mingw32-make -e CC=gcc
make install

cd ..\..\

cmake CMakeLists.txt
make -e CC=gcc