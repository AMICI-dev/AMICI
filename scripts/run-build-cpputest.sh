# Cpputest
mkdir -p ${AMICI_PATH}/ThirdParty
cd ${AMICI_PATH}/ThirdParty 
if [ ! -d "cpputest-3.8" ]; then
    if [ ! -e "cpputest-3.8.tar.gz" ]; then
        wget https://github.com/cpputest/cpputest/releases/download/v3.8/cpputest-3.8.tar.gz
    fi
    tar -xzf cpputest-3.8.tar.gz
    cd cpputest-3.8/cpputest_build/
    ../configure && make
fi
