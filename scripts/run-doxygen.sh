#!/bin/bash
# generate code documentation via doxygen

AMICI_PATH="`dirname \"$0\"`"
AMICI_PATH="`( cd \"$AMICI_PATH/..\" && pwd )`"

# read environment variables put there by build script
if [ -f ${AMICI_PATH}/scripts/env.sh ]; then
    . ${AMICI_PATH}/scripts/env.sh
fi

# build mtocpp
cd ${AMICI_PATH}/ThirdParty
if [ ! -d "mtocpp-master" ]; then
    if [ ! -e "mtocpp-master.zip" ]; then
        wget -O mtocpp-master.zip https://github.com/mdrohmann/mtocpp/archive/master.zip
    fi
    unzip mtocpp-master.zip
    mkdir ./mtocpp-master/build
    cd ./mtocpp-master/build && cmake .. && make
    if [ $? -ne 0 ] ; then
        exit 1
    fi
fi

cd ${AMICI_PATH}

# generate filter
echo "$AMICI_PATH/ThirdParty/mtocpp-master/build/mtocpp \$1 $AMICI_PATH/mtoc/config/mtocpp.conf" > ${AMICI_PATH}/mtoc/config/mtocpp_filter.sh

chmod +x ${AMICI_PATH}/mtoc/config/mtocpp_filter.sh

# generate doxyfile
cp ${AMICI_PATH}/mtoc/config/Doxyfile.template ${AMICI_PATH}/mtoc/config/Doxyfile

sed -i -e "s#_OutputDir_#$AMICI_PATH/doc#g" ${AMICI_PATH}/mtoc/config/Doxyfile
sed -i -e "s#_SourceDir_#$AMICI_PATH#g" ${AMICI_PATH}/mtoc/config/Doxyfile
sed -i -e "s#_ConfDir_#$AMICI_PATH/mtoc/config#g" ${AMICI_PATH}/mtoc/config/Doxyfile
sed -i -e "s#_ProjectName_#AMICI#g" ${AMICI_PATH}/mtoc/config/Doxyfile
sed -i -e "s#_ProjectDescription_#Advanced Matlab Interface for CVODES and IDAS#g" ${AMICI_PATH}/mtoc/config/Doxyfile
sed -i -e "s#_ProjectLogo_##g" ${AMICI_PATH}/mtoc/config/Doxyfile
sed -i -e "s#_ProjectVersion_##g" ${AMICI_PATH}/mtoc/config/Doxyfile
sed -i -e "s#_MTOCFILTER_#$AMICI_PATH/mtoc/config/mtocpp_filter.sh#g" ${AMICI_PATH}/mtoc/config/Doxyfile
sed -i -e "s#_LatexExtras_#$AMICI_PATH/mtoc/config/latexextras#g" ${AMICI_PATH}/mtoc/config/Doxyfile
sed -i -e "s#_GenLatex_#NO#g" ${AMICI_PATH}/mtoc/config/Doxyfile
sed -i -e "s#_HaveDot_#YES#g" ${AMICI_PATH}/mtoc/config/Doxyfile
sed -i -e "s#WARN_LOGFILE      =#WARN_LOGFILE      =$AMICI_PATH/mtoc/warnings.log#g" ${AMICI_PATH}/mtoc/config/Doxyfile

# generate latexextras

cp ${AMICI_PATH}/mtoc/config/latexextras.template ${AMICI_PATH}/mtoc/config/latexextras.sty
sed -i -e "s#_ConfDir_#$AMICI_PATH/mtoc/config#g" ${AMICI_PATH}/mtoc/config/latexextras.sty

doxygen "$AMICI_PATH/mtoc/config/Doxyfile"


#cleanup
#rm ${AMICI_PATH}/mtoc/config/latexextras.sty
rm ${AMICI_PATH}/mtoc/config/Doxyfile
rm ${AMICI_PATH}/mtoc/config/mtocpp_filter.sh

# check if warnings log was created
if [ -f ${AMICI_PATH}/mtoc/warnings.log  ]; then
    # check if warnings log is empty
if [ -s ${AMICI_PATH}/mtoc/warnings.log ]; then
        echo "DOXYGEN failed:"
        cat ${AMICI_PATH}/mtoc/warnings.log
        rm ${AMICI_PATH}/mtoc/warnings.log
        cat ${AMICI_PATH}/doc/_formulas.log
        exit 1
    else
        exit 0
    fi
else
    exit 1
fi

#$AMICI_PATH/ThirdParty/mtocpp-master/build/mtocpp_post "$AMICI_PATH/doc"

#cd ${AMICI_PATH}/doc/latex
#make
#mv ./refman.pdf ${AMICI_PATH}/AMICI_guide.pdf
