cd .\SuiteSparse\SuiteSparse_config

nmake library -e CC="gcc"

cd .\..\AMD

$CurrentDir = $(get-location).Path;

nmake library -e CC="gcc" 

cd .\..\BTF

nmake library -e CC="gcc" 

cd .\..\CAMD

nmake library -e CC="gcc" 

cd .\..\COLAMD

nmake library -e CC="gcc"  

cd .\..\KLU

nmake library -e CC="gcc" 


mkdir .\..\..\sundials\build\
cd .\..\..\sundials\build\

ls C:\projects\amici\SuiteSparse\lib

cmake .. -DCMAKE_INSTALL_PREFIX="C:/projects/amici/build/sundials" `
-DBUILD_ARKODE=OFF `
-DBUILD_CVODE=OFF `
-DBUILD_IDA=OFF `
-DBUILD_KINSOL=OFF `
-DBUILD_SHARED_LIBS=ON `
-DBUILD_STATIC_LIBS=OFF `
-DEXAMPLES_ENABLE=OFF `
-DEXAMPLES_INSTALL=OFF `
-DKLU_ENABLE=ON `
-DKLU_LIBRARY_DIR="C:/projects/amici/SuiteSparse/lib" `
-DKLU_INCLUDE_DIR="C:/projects/amici/SuiteSparse/include" `

msbuild ALL BUILD.vcxproj
msbuild INSTALL.vcxproj

cd ..\..

cmake CMakeLists.txt
msbuild model_dirac.vcxproj

ls

if ( (Test-Path ".\main.exe") -eq $false)
{
	throw "build unsuccessfull, model_dirac failed to compile!"
}