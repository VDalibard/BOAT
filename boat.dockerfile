FROM ubuntu:16.04

### get wget git etc
RUN apt-get update; apt-get -y install git
RUN apt-get update; apt-get -y install wget
RUN apt-get update; apt-get -y install gcc
RUN apt-get update; apt-get -y install g++
RUN apt-get update; apt-get -y install gfortran
RUN apt-get update; apt-get -y install cmake
RUN apt-get update; apt-get -y libblas-dev
RUN apt-get update; apt-get -y libblas-dev

# installing dependencies: boost
RUN apt-get update; apt-get -y install libboost-dev

# installing dependencies: eigen
WORKDIR "/var/tmp"
RUN git clone https://github.com/eigenteam/eigen-git-mirror
WORKDIR "/var/tmp/eigen-git-mirror"
RUN mkdir build
WORKDIR "/var/tmp/eigen-git-mirror/build"
RUN cmake ..; make; make install

# installing dependencies: nlopt
WORKDIR "/var/tmp"
RUN git clone git://github.com/stevengj/nlopt
WORKDIR "/var/tmp/nlopt"
RUN mkdir build
WORKDIR "/var/tmp/nlopt/build"
RUN cmake ..; make; make install

# installing BOAT
WORKDIR "/"
RUN git clone --branch chk_depend https://github.com/alan-turing-institute/BOAT.git

WORKDIR "/BOAT/examples/branin_hoo/"

CMD ["make"]
