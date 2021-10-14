# In this image we do not include the R-based plot scripts.
FROM ubuntu:20.04 as build

# Set noninterative mode
ENV DEBIAN_FRONTEND noninteractive
ENV LD_LIBRARY_PATH=/usr/local/lib/:$LD_LIBRARY_PATH

# apt update and install build requirements
RUN apt-get update \
     && apt-get install -y \
        g++=4:9.3.0-1ubuntu2 \
        cmake=3.16.3-1ubuntu1 \
        git \
        wget \
        libbz2-dev=1.0.8-2 \
        libcurl4-openssl-dev=7.68.0-1ubuntu2.7 \
        zlib1g-dev=1:1.2.11.dfsg-2ubuntu1.2 \
        liblzma-dev=5.2.4-1ubuntu1

# Compile htslib
WORKDIR /deps
RUN wget -q https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2 \
    && tar -xf htslib-1.11.tar.bz2 \
    && mv htslib-1.11 htslib
WORKDIR /deps/htslib
RUN autoheader; autoconf; ./configure --prefix=/usr/local/ \
    && make && make install

# Compile VerifyBamID. Version 2.0.1
WORKDIR /
RUN git clone --depth 1 --branch 2.0.1 git://github.com/Griffan/VerifyBamID.git
WORKDIR /VerifyBamID/build
RUN  cmake .. \
     && make \
     && make test

# Final image
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND noninteractive
RUN apt update && apt install -y \
        libgomp1=10.3.0-1ubuntu1~20.04 \
        libcurl4-openssl-dev=7.68.0-1ubuntu2.7 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /VerifyBamID/bin/VerifyBamID /usr/local/bin/
COPY --from=build /usr/local/lib/ /usr/local/lib/
CMD ["VerifyBamID"]
