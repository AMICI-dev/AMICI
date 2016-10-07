cd .\SuiteSparse\SuiteSparse_config

mingw32-make library -e CC="gcc"

cd .\..\AMD

mingw32-make library -e CC="gcc" 

cd .\..\BTF

mingw32-make library -e CC="gcc" 

cd .\..\CAMD

mingw32-make library -e CC="gcc" 

cd .\..\COLAMD

mingw32-make library -e CC="gcc"  

cd .\..\KLU

mingw32-make library -e CC="gcc" 

mkdir .\..\..\sundials\build\
cd .\..\..\sundials\build\

$envPaths = $env:Path -split ';'

if ($envPaths -contains 'C:\Program Files\Git\usr\bin') {
    $envPaths = $envPaths | where { $_ -and $_ -ne 'C:\Program Files\Git\usr\bin' }
    $env:Path = $envPaths -join ';'
}

cmake .. -DCMAKE_INSTALL_PREFIX="C:/projects/amici/build/sundials" `
-DBUILD_ARKODE=OFF `
-DBUILD_CVODE=OFF `
-DBUILD_IDA=OFF `
-DBUILD_KINSOL=OFF `
-DBUILD_SHARED_LIBS=OFF `
-DBUILD_STATIC_LIBS=ON `
-DEXAMPLES_ENABLE=OFF `
-DEXAMPLES_INSTALL=OFF `
-DKLU_ENABLE=ON `
-DKLU_LIBRARY_DIR="C:/projects/amici/SuiteSparse/lib" `
-DKLU_INCLUDE_DIR="C:/projects/amici/SuiteSparse/include" `
-G "MinGW Makefiles"

mingw32-make
mingw32-make install

cd ..\..

ls hdf5

Start-Process -FilePath msiexec -ArgumentList /i, "hdf5\HDF5-1.8.17-win32.msi", -Wait

ls hdf5

ls %APPVEYOR_BUILD_FOLDER%\libhdf5-lib

cmake CMakeLists.txt -G "MinGW Makefiles"

mingw32-make

ls

if ( (Test-Path ".\main.exe") -eq $false)
{
	throw "build unsuccessfull, model_dirac failed to compile!"
}