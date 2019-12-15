#!/bin/bash
SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd $SCRIPT_PATH/.. && pwd)

# avtivate virtual environment
source ${AMICI_PATH}/build/venv/bin/activate

# test PEtab command line import
MODEL_NAME="deleteme"
MODEL_OUTPUT_DIR="${MODEL_NAME}"
SBML_FILE_NAME="https://raw.githubusercontent.com/LeonardSchmiester/Benchmark-Models/hackathon/hackathon_contributions_new_data_format/Zheng_PNAS2012/model_Zheng_PNAS2012.xml"
MEASUREMENT_FILE_NAME="https://raw.githubusercontent.com/LeonardSchmiester/Benchmark-Models/hackathon/hackathon_contributions_new_data_format/Zheng_PNAS2012/measurementData_Zheng_PNAS2012.tsv"
CONDITION_FILE_NAME="https://raw.githubusercontent.com/LeonardSchmiester/Benchmark-Models/hackathon/hackathon_contributions_new_data_format/Zheng_PNAS2012/experimentalCondition_Zheng_PNAS2012.tsv"
PARAMETER_FILE_NAME="https://raw.githubusercontent.com/LeonardSchmiester/Benchmark-Models/hackathon/hackathon_contributions_new_data_format/Zheng_PNAS2012/parameters_Zheng_PNAS2012.tsv"
SBML_FILE_NAME2="deleteme_model.xml"
curl "${SBML_FILE_NAME}" > "${SBML_FILE_NAME2}"
amici_import_petab.py -o "${MODEL_OUTPUT_DIR}" --no-compile \
  -s "${SBML_FILE_NAME2}" -m "${MEASUREMENT_FILE_NAME}" \
  -c "${CONDITION_FILE_NAME}" -p "${PARAMETER_FILE_NAME}" -n "${MODEL_NAME}"
rm -rf "${MODEL_OUTPUT_DIR}"
