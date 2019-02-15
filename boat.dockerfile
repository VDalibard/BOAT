FROM ubuntu:16.04

### get wget git etc
RUN apt-get update; apt-get -y install git
RUN apt-get update; apt-get -y install wget
RUN apt-get update; apt-get -y install gcc g++ gfortran
RUN apt-get update; apt-get -y install cmake
RUN apt-get update; apt-get -y install libblas-dev
RUN apt-get update; apt-get -y install liblapack-dev
RUN apt-get update; apt-get -y install libboost-all-dev

# installing dependencies: eigen
WORKDIR "/var/tmp"
RUN git clone https://github.com/eigenteam/eigen-git-mirror
WORKDIR "/var/tmp/eigen-git-mirror"
RUN mkdir build
WORKDIR "/var/tmp/eigen-git-mirror/build"
RUN cmake ..; make; make install

#ENV PATH="/usr/local/include/eigen3:${PATH}"
ENV EIGEN_DIR="/usr/local/include/eigen3"
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"

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
RUN mkdir BOAT/build
WORKDIR "/BOAT/build/"
RUN cmake ../src; make

# WORKDIR "/BOAT/examples/branin_hoo/"
# CMD ["make"]

WORKDIR "/"
