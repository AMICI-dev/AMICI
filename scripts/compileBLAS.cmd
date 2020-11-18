echo compileBLAS.cmd started for openBLAS version %1
cd "C:\BLAS\OpenBLAS-v%1\OpenBLAS-%1"
call "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
cmake -G "Ninja" ^
    -DBUILD_DOUBLE=1 ^
    -DBUILD_SHARED_LIBS=ON ^
    -DCMAKE_INSTALL_PREFIX:PATH="C:\BLAS\OpenBLAS-v%1\OpenBLAS-%1\out\install\x64-Release" ^
    -DCMAKE_C_COMPILER:FILEPATH=cl ^
    -DCMAKE_BUILD_TYPE=Release ^
    -DCMAKE_MAKE_PROGRAM=ninja ^
    "C:\BLAS\OpenBLAS-v%1\OpenBLAS-%1"
cmake --build "C:\BLAS\OpenBLAS-v%1\OpenBLAS-%1" --parallel 2
echo compileBLAS.cmd completed
