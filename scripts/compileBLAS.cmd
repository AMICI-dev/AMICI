echo compileBLAS.cmd started for openBLAS version %1
cd /D "C:\BLAS\OpenBLAS-%1\OpenBLAS-%1"
call "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
cmake -S . -B build ^
    -G "Ninja" ^
    -DBUILD_DOUBLE=1 ^
    -DBUILD_SHARED_LIBS=ON ^
    -DCMAKE_INSTALL_PREFIX:PATH="C:/BLAS/OpenBLAS" ^
    -DCMAKE_C_COMPILER:FILEPATH=cl ^
    -DCMAKE_BUILD_TYPE=Release ^
    -DCMAKE_MAKE_PROGRAM=ninja
cmake --build build --parallel
cmake --install build
echo compileBLAS.cmd completed
