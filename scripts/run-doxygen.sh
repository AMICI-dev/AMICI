#!/bin/bash
# generate code documentation via doxygen

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

# build mtocpp
cd ${AMICI_PATH}/ThirdParty
if [ ! -d "mtocpp-master" ]; then
    if [ ! -e "mtocpp-master.zip" ]; then
        wget -O mtocpp-master.zip https://github.com/mdrohmann/mtocpp/archive/master.zip
    fi
    unzip mtocpp-master.zip
    mkdir ./mtocpp-master/build
    cd ./mtocpp-master/build && cmake .. && make mtocpp mtocpp_post
    if [ $? -ne 0 ] ; then
        exit 1
    fi
fi

cd ${AMICI_PATH}
MTOC_CONFIG_PATH=${AMICI_PATH}/matlab/mtoc/config
# generate filter
echo "$AMICI_PATH/ThirdParty/mtocpp-master/build/mtocpp \$1 ${MTOC_CONFIG_PATH}/mtocpp.conf" > ${MTOC_CONFIG_PATH}/mtocpp_filter.sh

chmod +x ${MTOC_CONFIG_PATH}/mtocpp_filter.sh

# generate doxyfile
DOXYFILE=${MTOC_CONFIG_PATH}/Doxyfile
cp ${MTOC_CONFIG_PATH}/Doxyfile.template ${DOXYFILE}
DOXY_WARNING_FILE=${AMICI_PATH}/matlab/mtoc/warnings.log

sed -i -e "s#_OutputDir_#$AMICI_PATH/doc#g" ${DOXYFILE}
sed -i -e "s#_SourceDir_#$AMICI_PATH#g" ${DOXYFILE}
sed -i -e "s#_ConfDir_#${MTOC_CONFIG_PATH}#g" ${DOXYFILE}
sed -i -e "s#_ProjectName_#AMICI#g" ${DOXYFILE}
sed -i -e "s#_ProjectDescription_#Advanced Multilanguage Interface for CVODES and IDAS#g" ${DOXYFILE}
sed -i -e "s#_ProjectLogo_##g" ${DOXYFILE}
sed -i -e "s#_ProjectVersion_##g" ${DOXYFILE}
sed -i -e "s#_MTOCFILTER_#${MTOC_CONFIG_PATH}/mtocpp_filter.sh#g" ${DOXYFILE}
sed -i -e "s#_LatexExtras_#${MTOC_CONFIG_PATH}/latexextras#g" ${DOXYFILE}
sed -i -e "s#_GenLatex_#YES#g" ${DOXYFILE}
sed -i -e "s#_HaveDot_#YES#g" ${DOXYFILE}
sed -i -e "s#WARN_LOGFILE      =#WARN_LOGFILE      =${DOXY_WARNING_FILE}#g" ${DOXYFILE}

# generate latexextras

cp ${MTOC_CONFIG_PATH}/latexextras.template ${MTOC_CONFIG_PATH}/latexextras.sty
sed -i -e "s#_ConfDir_#${MTOC_CONFIG_PATH}#g" ${MTOC_CONFIG_PATH}/latexextras.sty

export PATH=/Library/TeX/texbin:$PATH

doxygen "${DOXYFILE}"

#cleanup
#rm ${AMICI_PATH}/mtoc/config/latexextras.sty
rm ${DOXYFILE}
rm ${MTOC_CONFIG_PATH}/mtocpp_filter.sh

cd ${AMICI_PATH}/doc/latex

make

cat refman.log
cat refman.out

cp ./refman.pdf ${AMICI_PATH}/AMICI_guide.pdf

# suppress doxygen warnings about status badges
grep -v "warning: Unexpected html tag <img> found within <a href=...> context" ${DOXY_WARNING_FILE} > ${DOXY_WARNING_FILE}_tmp
mv ${DOXY_WARNING_FILE}_tmp ${DOXY_WARNING_FILE}

# suppress doxygen warning about unresolved external links (problem unclear)
grep -v "warning: unable to resolve reference to \`https" ${DOXY_WARNING_FILE} > ${DOXY_WARNING_FILE}_tmp
mv ${DOXY_WARNING_FILE}_tmp ${DOXY_WARNING_FILE}

grep -v "error: Problem running ghostscript gs -q -g146x60 -r384x384x -sDEVICE=ppmraw -sOutputFile=_form0.pnm -dNOPAUSE -dBATCH -- _form0.ps. Check your installation!" ${DOXY_WARNING_FILE} > ${DOXY_WARNING_FILE}_tmp
mv ${DOXY_WARNING_FILE}_tmp ${DOXY_WARNING_FILE}

# check if warnings log was created
if [ -f ${DOXY_WARNING_FILE}  ]; then
    # check if warnings log is empty
    if [ -s ${DOXY_WARNING_FILE} ]; then
        echo "DOXYGEN failed:"
        cat ${DOXY_WARNING_FILE}
        #rm ${DOXY_WARNING_FILE}
        exit 1
    else
        exit 0
    fi
else
    exit 1
fi
