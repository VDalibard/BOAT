FROM ubuntu:16.04

### get wget git etc
RUN apt-get update; apt-get -y install git
RUN apt-get update; apt-get -y install wget
RUN apt-get update; apt-get -y install gcc
RUN apt-get update; apt-get -y install cmake

# installing dependencies
RUN apt-get update; apt-get -y install libboost-dev
RUN apt-get update; apt-get -y install libeigen3-dev

# installing dependencies: nlopt
WORKDIR "/var/tmp"
RUN git clone git://github.com/stevengj/nlopt
WORKDIR "/var/tmp/nlopt"
RUN mkdir build
WORKDIR "/var/tmp/nlopt/build"
RUN cmake ..; make; make install

# installing BOAT
WORKDIR "/"
RUN git clone https://github.com/alan-turing-institute/BOAT.git
