#!/bin/bash
# Download and build mtocpp (Doxygen filter for Matlab)

set -euo pipefail
set -x

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd "$SCRIPT_PATH/.." && pwd)

MTOC_CONFIG_PATH=${AMICI_PATH}/matlab/mtoc/config

cd "${AMICI_PATH}/ThirdParty"

# download mtocpp?
if [ ! -d "mtocpp-master" ]; then
    # mtocpp requires ragel
    if ! command -v ragel &> /dev/null; then
      echo "ragel not found"

      if [ ! -d "ragel-6.10" ]; then
        if [ ! -e "ragel-6.10.tar.gz" ]; then
            echo "Downloading ragel ..."
            wget -O ragel-6.10.tar.gz http://www.colm.net/files/ragel/ragel-6.10.tar.gz
        fi
        # build ragel?
        tar -xzf ragel-6.10.tar.gz
        (cd ragel-6.10 && ./configure --prefix="$(pwd)/install" && make -j2 && make install)
      fi
      export PATH=$PATH:$(pwd)/ragel-6.10/install/bin
    fi

    if [ ! -e "mtocpp-master.zip" ]; then
        echo "Downloading mtocpp ..."
        wget -O mtocpp-master.zip https://github.com/mdrohmann/mtocpp/archive/master.zip
    fi
    # build mtocpp?
    unzip mtocpp-master.zip

    # patch for xml support for postprocessor
    sed -i.bak 's/== "tex"$/== "tex" || file.substr(file.find_last_of(".") + 1) == "xml"/' mtocpp-master/src/postprocess.rl

    mkdir -p mtocpp-master/build

    if command -v cmake &> /dev/null; then
      echo "Building mtocpp using CMake..."
      cd mtocpp-master/build && cmake .. && make mtocpp mtocpp_post
      if [ $? -ne 0 ] ; then
          exit 1
      fi
    else
      # No CMake on ReadTheDocs :(
      echo "Building mtocpp without CMake..."
      cd mtocpp-master/src
      sed 's/@MTOC++_VERSION_MAJOR@/1/;s/@MTOC++_VERSION_MINOR@/5/' < config.h.in > config.h
      ragel -C -T0 -o confscanner.cc confscanner.rl
      ragel -C -T0 -o mfilescanner_parser.cc mfilescanner_parser.rl
      c++   -I"$(pwd)" -o mtocpp.cc.o -c mtocpp.cc
      c++   -I"$(pwd)" -o mfilescanner_parser.cc.o -c mfilescanner_parser.cc
      c++   -I"$(pwd)" -o mfilescanner.cc.o -c mfilescanner.cc
      c++   -I"$(pwd)" -o confscanner.cc.o -c confscanner.cc
      c++   -rdynamic mtocpp.cc.o mfilescanner_parser.cc.o mfilescanner.cc.o confscanner.cc.o -o mtocpp
      ragel -C -T0 -o postprocess.cc postprocess.rl
      c++   -I"$(pwd)" -o postprocess.cc.o -c postprocess.cc
      c++   -rdynamic postprocess.cc.o -o mtocpp_post
      cd ../build/
      ln -s ../src/mtocpp mtocpp
      ln -s ../src/mtocpp_post mtocpp_post
    fi
fi

# generate filter
echo "$AMICI_PATH/ThirdParty/mtocpp-master/build/mtocpp \$1 ${MTOC_CONFIG_PATH}/mtocpp.conf" > "${MTOC_CONFIG_PATH}/mtocpp_filter.sh"
chmod +x "${MTOC_CONFIG_PATH}/mtocpp_filter.sh"
