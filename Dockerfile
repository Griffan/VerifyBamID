# In this image we do not include the R-based plot scripts.
FROM ubuntu:26.04 AS build

# Set noninterative mode
ENV DEBIAN_FRONTEND=noninteractive

# The release to build, and the repository to clone it from. VERSION must name
# a release whose build bundles htslib.
ARG VERSION=v2.0.3
ARG REPO=https://github.com/Griffan/VerifyBamID.git

# apt update and install build requirements. The build downloads htslib and
# libdeflate over HTTPS and links them statically.
RUN apt-get update \
     && apt-get install -y \
        build-essential \
        ca-certificates \
        cmake \
        g++ \
        git \
        libbz2-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        zlib1g-dev \
        liblzma-dev

# Compile VerifyBamID. Version ${VERSION}
WORKDIR /
RUN git clone --depth 1 --branch ${VERSION} ${REPO} VerifyBamID
WORKDIR /VerifyBamID/build

# -fsigned-char: VerifyBamID's genotype-missing sentinel relies on `char`
# being signed (default on x86_64 but not on aarch64), which otherwise
# breaks missing-genotype detection and fails testReadVcf.
RUN cmake -DCMAKE_C_FLAGS=-fsigned-char -DCMAKE_CXX_FLAGS=-fsigned-char .. && \
    make -j"$(nproc)" && \
    make test

# Final image
FROM ubuntu:26.04
ENV DEBIAN_FRONTEND="noninteractive"

RUN apt update && apt install -y \
        ca-certificates \
        libgomp1 \
        libcurl4t64 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /VerifyBamID/bin/VerifyBamID /usr/local/bin/

CMD ["VerifyBamID"]
