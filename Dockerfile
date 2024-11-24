FROM ubuntu:20.04

RUN apt-get -y update
RUN apt-get -y install less git

COPY jdk-6u45-linux-x64.bin /jdk-6u45-linux-x64.bin
RUN bash /jdk-6u45-linux-x64.bin

RUN git clone https://github.com/smirarab/sepp.git

RUN apt-get -y install wget
RUN wget https://archive.apache.org/dist/ant/binaries/apache-ant-1.7.1-bin.tar.gz
RUN tar xzvf apache-ant-1.7.1-bin.tar.gz 

ENV PATH=/apache-ant-1.7.1/bin:/jdk1.6.0_45/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
