FROM ubuntu:24.04

# install dependencies

RUN apt-get update && apt-get install -y \
    g++ \
    git \
    libboost-serialization-dev \
    libboost-chrono-dev \
    libhdf5-serial-dev \
    python-is-python3 \
    python3 \
    python3-dev \
    python3-pip \
    python3-venv \
    wget

# BEGIN BINDER
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER=${NB_USER}
ENV NB_UID=${NB_UID}
ENV HOME=/home/${NB_USER}
ENV VENV_PATH=${HOME}/venv
ENV PATH="${VENV_PATH}/bin:${PATH}"
# If there is no default 1000 user:
#RUN adduser --disabled-password \
#    --gecos "Default user" \
#    --uid ${NB_UID} \
#    ${NB_USER}
# On ubuntu>=24.04 images there is ubuntu(1000)
#  see https://bugs.launchpad.net/cloud-images/+bug/2005129
RUN usermod -l ${NB_USER} ubuntu && \
    mv /home/ubuntu ${HOME} && \
    usermod -d ${HOME} ${NB_USER} && \
    chown -R ${NB_UID}:${NB_UID} ${HOME}
USER ${NB_USER}
RUN python3 -m venv ${VENV_PATH} && \
    ${VENV_PATH}/bin/pip install --upgrade build wheel setuptools pip && \
    ${VENV_PATH}/bin/pip install --no-cache-dir notebook jupyterlab jupyterhub jax[cpu] antimony
# END BINDER

# BEGIN BINDER
# Make sure the contents of our repo are in ${HOME}
USER root
COPY . ${HOME}
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}
# END BINDER

WORKDIR ${HOME}
RUN . ./.profile && ${VENV_PATH}/bin/python3 -m build --sdist python/sdist && \
    ${VENV_PATH}/bin/pip install git+https://github.com/pysb/pysb@master && \
    ${VENV_PATH}/bin/pip install -U sympy && \
    ${VENV_PATH}/bin/pip install -v $(ls -t python/sdist/dist/amici-*.tar.gz | head -1)[petab,vis,jax,sciml] && \
    ${VENV_PATH}/bin/pip install "git+https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab.git@master#subdirectory=src/python&egg=benchmark_models_petab" && \
    scripts/buildBNGL.sh

ENV BNGPATH="${HOME}/ThirdParty/BioNetGen-2.7.0"
