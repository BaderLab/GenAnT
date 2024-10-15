#docker file for docker to run mikado 2.3.2 - modified from docker file from - https://github.com/EI-CoreBioinformatics/mikado/blob/master/Docker/Dockerfile.ubuntu
FROM ubuntu:22.04
ENV PYTHONDONTWRITEBYTECODE true

#samtools version
ENV hts_version 1.11

#portucullis version
ENV port_version 1.2.2

#diamond linux version
ENV dmd_version 2.0.6

#Prodigal version
ENV prod_version 2.6.3

WORKDIR /usr/local/src

#install required software
RUN apt update && apt install -y software-properties-common && \
    add-apt-repository universe && add-apt-repository multiverse 
RUN add-apt-repository -y restricted && apt update
RUN apt install -y zlib1g-dev libboost-dev '^libboost-.*71-dev' build-essential automake autoconf curl 

#install required python and packages
RUN apt install -y libbz2-dev libncurses-dev lzma-dev python3 python3-numpy python3-dev python3-scipy 
RUN DEBIAN_FRONTEND=noninteractive apt-get -y install python3-pandas python3-pip liblzma-dev libtool
RUN pip3 install wheel && pip3 install sqlalchemy==1.4.41 && pip3 install sqlalchemy-utils==0.38.3 && pip3 install python-rapidjson && pip3 install "Mikado==2.3.2" 

#install samtools
WORKDIR /usr/local/src
RUN curl -L https://github.com/samtools/samtools/releases/download/${hts_version}/samtools-${hts_version}.tar.bz2 | tar xj
RUN cd samtools-${hts_version} && autoreconf && ./configure && make -j 10 && make install 
RUN cd /usr/local/src 

#install portcullis
RUN curl -L https://github.com/maplesond/portcullis/archive/refs/tags/${port_version}.tar.gz | tar -xz 
#RUN cd /usr/local/src/portcullis-${port_version} && ./autogen.sh && ./configure && make -j 10 && make -j 10 check && make install 
#RUN mkdir -p /usr/local/src/ 
RUN cd /usr/local/src 
RUN curl -L https://github.com/bbuchfink/diamond/releases/download/v${dmd_version}/diamond-linux64.tar.gz | tar -xz -C /usr/local/bin/ 

#install Prodigal
RUN curl -L https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux > /usr/local/bin/prodigal 
RUN chmod +x /usr/local/bin/prodigal 
RUN cd /usr/local/ 
RUN rm -rf src/*

#install sqlite3 
RUN apt install sqlite3

WORKDIR /global
CMD ["/bin/bash"]
