echo compileBLAS.cmd started
cd "C:\BLAS\OpenBLAS-v0.3.10\OpenBLAS-0.3.10"
call "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
cmake -G "Ninja" -DBUILD_DOUBLE -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX:PATH="C:\BLAS\OpenBLAS-v0.3.10\OpenBLAS-0.3.10\out\install\x64-Release" -DCMAKE_C_COMPILER:FILEPATH=cl -DCMAKE_BUILD_TYPE=Release -DCMAKE_MAKE_PROGRAM=ninja "C:\BLAS\OpenBLAS-v0.3.10\OpenBLAS-0.3.10"
cmake --build --parallel "C:\BLAS\OpenBLAS-v0.3.10\OpenBLAS-0.3.10"
echo compileBLAS.cmd completed
