# Source Image
 FROM ubuntu:latest
 FROM gcc:7.5
  # Set noninterative mode
 ENV DEBIAN_FRONTEND noninteractive

  # apt update and install global requirements
 RUN apt-get clean all && \
     apt-get update && \
     apt-get install -y  \
         autoconf \
         cmake \
         git \
         libbz2-dev \
         libcurl4-openssl-dev \
         libssl-dev \
         zlib1g-dev \
         liblzma-dev

  # apt clean and remove cached source lists
 RUN apt-get clean && \
     rm -rf /var/lib/apt/lists/*

 RUN git clone git://github.com/samtools/htslib.git
 RUN cd htslib && \
     autoheader && \
     autoconf && \
     ./configure --prefix=/usr/local/ && \
     make && \
     make install

 RUN git clone git://github.com/Griffan/VerifyBamID.git
 RUN cd VerifyBamID && \
     mkdir build && \
     cd build && \
     cmake .. && \
     make && \
     make test
 RUN cp /VerifyBamID/bin/VerifyBamID /usr/local/bin

  # Define default command
 CMD ["VerifyBamID"]
