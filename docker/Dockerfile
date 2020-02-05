FROM alpine

RUN apk add blas-dev \
    g++ \
    musl-dev \
    py3-numpy \
    python3 \
    python3-dev \
    swig

COPY amici.tar.gz /tmp

RUN pip3 install -U --upgrade pip \
    && mkdir -p /tmp/amici/ \
    && cd /tmp/amici \
    && tar -xzf ../amici.tar.gz \
    && cd /tmp/amici/python/sdist \
    && python3 setup.py -v sdist \
    && pip3 install -v $(ls -t /tmp/amici/python/sdist/dist/amici-*.tar.gz | head -1)[petab] \
    && rm -rf /tmp && mkdir /tmp
