# Source Image
 FROM ubuntu:latest
# Set noninterative mode
 ENV DEBIAN_FRONTEND noninteractive
 ENV LD_LIBRARY_PATH=/usr/local/lib/:$LD_LIBRARY_PATH

  # apt update and install global requirements
 RUN apt-get clean all && \
     apt-get update && \
     apt-get install -y  \
         gcc-7 \
         g++-7 \
         autoconf \
         cmake \
         git \
         libbz2-dev \
         libcurl4-openssl-dev \
         zlib1g-dev \
         liblzma-dev

 RUN ln -sf /usr/bin/g++-7 /usr/bin/g++ && \
     ln -sf /usr/bin/gcc-7 /usr/bin/gcc
  # apt clean and remove cached source lists
 RUN apt-get clean && \
     rm -rf /var/lib/apt/lists/*

 RUN git clone git://github.com/samtools/htslib.git && \
     cd htslib && \
     autoheader && \
     autoconf && \
     ./configure --prefix=/usr/local/ && \
     make && \
     make install

 RUN git clone git://github.com/Griffan/VerifyBamID.git && \
     cd VerifyBamID && \
     mkdir build && \
     cd build && \
     cmake .. && \
     make && \
     make test

 RUN cp /VerifyBamID/bin/VerifyBamID /usr/local/bin

  # Define default command
 CMD ["VerifyBamID"]
