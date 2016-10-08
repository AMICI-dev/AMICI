#!/bin/bash
#cd 'C:/Program Files/HDF_Group/HDF5/1.8.17/bin/'
#for dll in *.dll; do
#	  base_def=$(basename $dll .dll)
#      def_file="${base_def}.def"
#      lib_file="lib${base_def}.a"
#      pexports $dll > $def_file
#      dlltool -d $def_file -l $lib_file
#done 
patch -p1 -d 'C:/Program Files/HDF_Group/HDF5/1.8.17/' < C:/projects/amici/scripts/hdf5-1.8.7-mingw.patch