#!/bin/bash
# Create a virtual environment and perform an editable amici installation
set -e

SCRIPT_PATH=$(dirname $BASH_SOURCE)
AMICI_PATH=$(cd "$SCRIPT_PATH/.." && pwd)

# Disabled until cmake package is made compatible with updated setup.py
#make python-wheel
#pip3 install --user --prefix= `ls -t ${AMICI_PATH}/build/python/amici-*.whl | head -1`

# test install from setup.py
venv_dir="${AMICI_PATH}/venv"
set +e
mkdir -p "${venv_dir}"
python3 -m venv "${venv_dir}" --clear
# in case this fails (usually due to missing ensurepip, try getting pip
# manually
if [[ $? ]]; then
    set -e
    python3 -m venv "${venv_dir}" --clear --without-pip
    source "${venv_dir}/bin/activate"
    get_pip=${AMICI_PATH}/get-pip.py
    curl "https://bootstrap.pypa.io/get-pip.py" -o "${get_pip}"
    python3 "${get_pip}"
    rm "${get_pip}"
else
    set -e
    source "${venv_dir}/bin/activate"
fi

# set python executable for cmake
export PYTHON_EXECUTABLE="${AMICI_PATH}/venv/bin/python"

python -m pip install --upgrade pip
DEP_GROUPS="${AMICI_PATH}/python/sdist/pyproject.toml"
python -m pip install \
  --group "${DEP_GROUPS}:build" \
  --group "${DEP_GROUPS}:pysb-dev" \
  --group "${DEP_GROUPS}:diffrax-dev" \
  --group "${DEP_GROUPS}:libpetab" \
  --group "${DEP_GROUPS}:petab-sciml"
python -m pip install optax # for jax petab notebook
AMICI_BUILD_TEMP="${AMICI_PATH}/python/sdist/build/temp" \
  python -m pip install --verbose \
    --group "${DEP_GROUPS}:test" \
    -e "${AMICI_PATH}/python/sdist[petab,vis,jax]" --no-build-isolation
python -m pip install --group "${DEP_GROUPS}:fiddy"
deactivate
