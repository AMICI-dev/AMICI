FROM ubuntu:24.04

# install dependencies

RUN apt-get update && apt-get install -y \
    g++ \
    libatlas-base-dev \
    libboost-serialization-dev \
    libboost-chrono-dev \
    libhdf5-serial-dev \
    python-is-python3 \
    python3 \
    python3-dev \
    python3-pip \
    python3-venv \
    swig

# prepare python install

COPY amici.tar.gz /from_host/
RUN bash -c 'python3 -m venv venv && \
    echo "source $HOME/venv/bin/activate" >> ~/.bashrc && \
    source venv/bin/activate && \
    pip3 install --upgrade pip build && \
    mkdir -p /from_host/amici/ && \
    cd /from_host/amici && \
    tar -xzf ../amici.tar.gz && cd python/sdist && \
    python3 -m build --sdist && \
    pip3 install -v $(ls -t dist/amici-*.tar.gz | head -1)[petab] && \
    rm -rf /from_host'
