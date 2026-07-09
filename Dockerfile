# In this image we do not include the R-based plot scripts.
FROM ubuntu:20.04 AS build

# Set noninterative mode
ENV DEBIAN_FRONTEND=noninteractive
ENV LD_LIBRARY_PATH="/usr/local/lib"

ARG HTSLIB_VERSION=1.23.1
ARG VERSION=v2.0.3

# apt update and install build requirements
RUN apt-get update \
     && apt-get install -y \
        build-essential \
        cmake \
        g++ \
        git \
        wget \
        libbz2-dev \
        libcurl4-openssl-dev \
        zlib1g-dev \
        liblzma-dev

# Compile htslib
WORKDIR /deps
RUN wget -q https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
    && tar -xf htslib-${HTSLIB_VERSION}.tar.bz2 \
    && mv htslib-${HTSLIB_VERSION} htslib

WORKDIR /deps/htslib
RUN autoheader; autoconf; ./configure --prefix=/usr/local/ \
    && make && make install

# Compile VerifyBamID. Version ${VERSION}
WORKDIR /
RUN git clone --depth 1 --branch ${VERSION} https://github.com/Griffan/VerifyBamID.git
WORKDIR /VerifyBamID/build

# -fsigned-char: VerifyBamID's genotype-missing sentinel relies on `char`
# being signed (default on x86_64 but not on aarch64), which otherwise
# breaks missing-genotype detection and fails testReadVcf.
RUN cmake -DCMAKE_C_FLAGS=-fsigned-char -DCMAKE_CXX_FLAGS=-fsigned-char .. && \
    make && \
    make test

# Final image
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND="noninteractive"

RUN apt update && apt install -y \
        libgomp1 \
        libcurl4-openssl-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /VerifyBamID/bin/VerifyBamID /usr/local/bin/
COPY --from=build /usr/local/lib/ /usr/local/lib/

CMD ["VerifyBamID"]
