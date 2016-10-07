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

mingw-get install mingw-utils

cd "C:/Program Files (x86)/HDF_Group/HDF5/1.8.17/lib"

for dll in `ls dll/*dll`; do
  def_file=`basename $dll .dll`.def
  lib_file=lib`basename $dll dll.dll`.a
  pexports $dll > $def_file
  dlltool -d $def_file -l lib/$lib_file
done

cd "C:/Program Files (x86)/HDF_Group/HDF5/1.8.17/include"

patch -p1 -d /c/swarm < C:/projects/amici/scripts/hdf5-1.8.7-mingw.patch

cmake CMakeLists.txt -G "MinGW Makefiles"

mingw32-make 

ls

if ( (Test-Path ".\main.exe") -eq $false)
{
	throw "build unsuccessfull, model_dirac failed to compile!"
}