$envPaths = $env:Path -split ';'

if ($envPaths -contains 'C:\MinGW\bin') {
    $envPaths = $envPaths | where { $_ -and $_ -ne 'C:\MinGW\bin' }
    $env:Path = $envPaths -join ';'
}

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

$envPaths = $env:Path -split ';'

if ($envPaths -contains 'C:\MinGW\msys\1.0\bin') {
    $envPaths = $envPaths | where { $_ -and $_ -ne 'C:\MinGW\msys\1.0\bin' }
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

ls "C:/Program Files/HDF_Group/HDF5/1.8.17/bin/"
ls "C:/Program Files/HDF_Group/HDF5/1.8.17/"
ls "C:/Program Files/HDF_Group/HDF5/1.8.17/cmake"
ls "C:/Program Files/HDF_Group/HDF5/1.8.17/lib"
ls "C:/Program Files/HDF_Group/HDF5/1.8.17/include"

cmake CMakeLists.txt `
-DCMAKE_VERBOSE_MAKEFILE=ON `
-G "MinGW Makefiles"

mingw32-make 

if ( (Test-Path ".\main.exe") -eq $false)
{
	throw "build unsuccessfull, model_dirac failed to compile!"
}