FROM ubuntu:20.04

RUN apt-get -y update
RUN apt-get -y install wget git

# you need to register an account with Oracle to legally obtain this file from https://www.oracle.com/java/technologies/javase-java-archive-javase6-downloads.html
COPY jdk-6u45-linux-x64.bin /jdk-6u45-linux-x64.bin
RUN bash /jdk-6u45-linux-x64.bin

RUN git clone https://github.com/smirarab/sepp.git

# download and "install" ant in version 1.7.1
RUN wget https://archive.apache.org/dist/ant/binaries/apache-ant-1.7.1-bin.tar.gz
RUN tar xzvf apache-ant-1.7.1-bin.tar.gz 

# ensure javac and ant can be found
ENV PATH=/apache-ant-1.7.1/bin:/jdk1.6.0_45/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# automatically compile and package seppJsonMerger.jar
CMD cd /sepp/tools/merge && ant compile && ant jar
