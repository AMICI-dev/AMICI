#!/bin/bash
# Generate code documentation via doxygen
set -euo pipefail
set -x

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

OUTDIR="${AMICI_PATH}/doc/build_doxygen"
mkdir -p "$OUTDIR"

cd "${AMICI_PATH}"

# Generate doxyfile
DOXYFILE="${AMICI_PATH}/doc/Doxyfile"
cp "${AMICI_PATH}/doc/Doxyfile.template" "${DOXYFILE}"
DOXY_WARNING_FILE=${AMICI_PATH}/doc/warnings.log

# Replace some template values
sed -i -e "s#_OutputDir_#$OUTDIR#g" ${DOXYFILE}
sed -i -e "s#_SourceDir_#$AMICI_PATH#g" ${DOXYFILE}
sed -i -e "s#_ConfDir_#${AMICI_PATH}/doc#g" ${DOXYFILE}
sed -i -e "s#_ProjectName_#AMICI#g" ${DOXYFILE}
sed -i -e "s#_ProjectDescription_#Advanced Multilanguage Interface for CVODES and IDAS#g" ${DOXYFILE}
sed -i -e "s#_ProjectLogo_##g" ${DOXYFILE}
sed -i -e "s#_ProjectVersion_##g" ${DOXYFILE}
sed -i -e "s#_GenLatex_#NO#g" ${DOXYFILE}
sed -i -e "s#_HaveDot_#YES#g" ${DOXYFILE}
# Fail if no replace was made
sed -i -re "/WARN_LOGFILE(\s*)=.*/{s##WARN_LOGFILE\1= ${DOXY_WARNING_FILE}#g;h};\${x;/./{x;q0};x;q1}" ${DOXYFILE}

# Run doxygen
doxygen "${DOXYFILE}"

#cleanup
rm "${DOXYFILE}"

# Build pdf
#cd "${OUTDIR}/latex"
#make
#cp ./refman.pdf "${AMICI_PATH}/AMICI_guide.pdf"

# suppress doxygen warnings about status badges
grep -v "warning: Unexpected html tag <img> found within <a href=...> context" "${DOXY_WARNING_FILE}" > "${DOXY_WARNING_FILE}_tmp" || [[ $? == 1 ]]
mv "${DOXY_WARNING_FILE}_tmp" "${DOXY_WARNING_FILE}"

# suppress doxygen warning about unresolved external links (problem unclear)
grep -v "warning: unable to resolve reference to \`https" "${DOXY_WARNING_FILE}" > "${DOXY_WARNING_FILE}_tmp" || [[ $? == 1 ]]
mv "${DOXY_WARNING_FILE}_tmp" "${DOXY_WARNING_FILE}"

grep -v "error: Problem.* running g.*. Check your installation!" "${DOXY_WARNING_FILE}" > "${DOXY_WARNING_FILE}_tmp" || [[ $? == 1 ]]
mv "${DOXY_WARNING_FILE}_tmp" "${DOXY_WARNING_FILE}"

# unclear error parsing ghostscript-generated eps files
# result seems fine despite error
grep -v "Couldn't extract bounding box from" "${DOXY_WARNING_FILE}" > "${DOXY_WARNING_FILE}_tmp" || [[ $? == 1 ]]
mv "${DOXY_WARNING_FILE}_tmp" "${DOXY_WARNING_FILE}"

# check if warnings log was created
if [ -f "${DOXY_WARNING_FILE}"  ]; then
    # check if warnings log is empty
    if [ -s "${DOXY_WARNING_FILE}" ]; then
        echo "DOXYGEN failed:"
        cat "${DOXY_WARNING_FILE}"
        #rm ${DOXY_WARNING_FILE}
        exit 1
    else
        exit 0
    fi
else
    exit 1
fi
