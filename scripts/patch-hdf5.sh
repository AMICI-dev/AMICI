#!/bin/bash
cd C:/Program Files (x86)/HDF_Group/HDF5/1.8.17/bin/
for dll in 'ls dll/*dll'; do
      def_file='basename $dll .dll'.def
      lib_file=lib'basename $dll dll.dll'.a
      pexports $dll > $def_file
      dlltool -d $def_file -l lib/$lib_file
done 
patch -p1 -d 'C:/Program Files (x86)/HDF_Group/HDF5/1.8.17/' < C:/projects/amici/scripts/hdf5-1.8.7-mingw.patch