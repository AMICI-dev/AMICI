This directory contains include files that are required for the Matlab-based installation without CMake.
Normally, those files would be generated by CMake. We can't add them to the other includes, because they 
would be found before the CMake-generated ones. Then, they will either not work on some systems, or they 
would disable some features we'd want.
