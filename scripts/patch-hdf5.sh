#!/bin/bash
cd 'C:/Program Files (x86)/HDF_Group/HDF5/1.8.17/bin/'
for dll in *.dll; do
	  echo "$dll"
	  base_def=$(basename $dll .dll)
	  echo "$base_def"
      def_file="${base_def}.def"
      base_lib=$(basename $dll dll.dll)
      echo "$base_lib"
      lib_file="lib${base_lib}.a"
      pexports $dll > $def_file
      dlltool -d $def_file -l lib/$lib_file
done 
patch -p1 -d 'C:/Program Files (x86)/HDF_Group/HDF5/1.8.17/' < C:/projects/amici/scripts/hdf5-1.8.7-mingw.patch