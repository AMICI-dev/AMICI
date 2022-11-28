FROM ubuntu:22.04

# install dependencies 

RUN apt update
RUN apt-get install -y \
    g++ \
    libatlas-base-dev \
    python3 \
    python3-dev \
    python3-pip \
    python3-venv \
    swig \
    # required for wheel-less libsbml installation 
    cmake

# prepare python install

COPY amici.tar.gz /tmp

RUN pip3 install --upgrade pip build

# build sdist

RUN mkdir -p /tmp/amici/ && \
    cd /tmp/amici && \
    tar -xzf ../amici.tar.gz && cd /tmp/amici/python/sdist && \
    python3 -m build --sdist

# install

RUN pip3 install -v $(ls -t /tmp/amici/python/sdist/dist/amici-*.tar.gz | head -1)[petab]

# cleanup

RUN rm -rf /tmp/*
