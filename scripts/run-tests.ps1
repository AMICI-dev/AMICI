cd .\SuiteSparse\SuiteSparse_config

mingw32-make library -e CC="gcc" LDFLAGS="$LDFLAGS -shared" 

cd .\..\AMD

mingw32-make library -e CC="gcc" LDFLAGS="$LDFLAGS -shared -L$SUITESPARSE/lib" 

cd .\..\BTF

mingw32-make library -e CC="gcc" LDFLAGS="$LDFLAGS -shared -L$SUITESPARSE/lib" 

cd .\..\CAMD

mingw32-make library -e CC="gcc" LDFLAGS="$LDFLAGS -shared -L$SUITESPARSE/lib" 

cd .\..\COLAMD

mingw32-make library -e CC="gcc" LDFLAGS="$LDFLAGS -shared -L$SUITESPARSE/lib" 

cd .\..\KLU

mingw32-make library -e CC="gcc" LDFLAGS="$LDFLAGS -shared -L$SUITESPARSE/lib" 


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
.\.. 

mingw32-make -e CC="gcc" LDFLAGS="$LDFLAGS -shared" 
make install

cd ..\..\

cmake CMakeLists.txt
mingw32-make -e CC=gcc