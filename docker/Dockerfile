FROM alpine

RUN apk add --update-cache --no-cache \
    blas-dev \
    g++ \
    libgcc \
    libquadmath \
    musl-dev \
    py3-numpy \
    python3 \
    python3-dev \
    swig
RUN echo "http://dl-8.alpinelinux.org/alpine/edge/testing" >> /etc/apk/repositories && apk add py3-pandas py3-matplotlib py3-scipy py3-h5py
#RUN addgroup user && useradd --no-log-init -g user user
#USER user
COPY amici.tar.gz /tmp

RUN pip3 install -U --upgrade pip wheel \
    && mkdir -p /tmp/amici/ \
    && cd /tmp/amici \
    && tar -xzf ../amici.tar.gz \
    && cd /tmp/amici/python/sdist \
    && python3 setup.py -v sdist \
    && pip3 install -v $(ls -t /tmp/amici/python/sdist/dist/amici-*.tar.gz | head -1)[petab] \
    && rm -rf /tmp && mkdir /tmp
