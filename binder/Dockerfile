FROM ubuntu:24.04

# install dependencies

RUN apt-get update && apt-get install -y \
    g++ \
    git \
    libatlas-base-dev \
    libboost-serialization-dev \
    libboost-chrono-dev \
    libhdf5-serial-dev \
    python-is-python3 \
    python3 \
    python3-dev \
    python3-pip \
    python3-venv \
    swig \
    wget

# BEGIN BINDER
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}
ENV PATH="${PATH}:$HOME/.local/bin"
RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
USER ${NB_USER}
RUN python3 -m pip install --upgrade build wheel setuptools pip && \
    python3 -m pip install --no-cache-dir notebook jupyterlab jupyterhub jax[cpu]
# END BINDER

# BEGIN BINDER
# Make sure the contents of our repo are in ${HOME}
USER root
COPY . ${HOME}
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}
# END BINDER

WORKDIR ${HOME}
RUN . ./.profile && python3 -m build --sdist python/sdist && \
    python3 -m pip install -v $(ls -t python/sdist/dist/amici-*.tar.gz | head -1)[petab,pysb,vis] && \
    python3 -m pip install "git+https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab.git@master#subdirectory=src/python&egg=benchmark_models_petab" && \
    scripts/buildBNGL.sh

ENV BNGPATH="${HOME}/ThirdParty/BioNetGen-2.7.0"
