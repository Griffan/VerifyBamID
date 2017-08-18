# Source Image
FROM ubuntu:latest

# Set noninterative mode
ENV DEBIAN_FRONTEND noninteractive

# apt update and install global requirements
RUN apt-get clean all && \
    apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y  \
        autoconf \
        build-essential \
        cmake \
        git \
        libbz2-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        zlib1g-dev

# apt clean and remove cached source lists
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install VerifyBamID
RUN git clone https://github.com/Griffan/VerifyBamID.git
RUN cd VerifyBamID && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make test
RUN cp /VerifyBamID/bin/VerifyBamID /usr/local/bin

# Define default command
CMD ["VerifyBamID"]
