CC=g++-8
CFLAGS=-Wall -fexceptions -std=c++14 -O3 -fopenmp
BOAT_DIR=./../../
LIBS= -lboat -lpthread -lblas -lnlopt -lm -lboost_serialization -lboost_system -lboost_coroutine -I $(EIGEN_DIR)

bo_naive_opti: bo_naive_opti.cpp
	$(CC) -o bo_naive_opti $(CFLAGS) -I$(BOAT_DIR)/include/ bo_naive_opti.cpp -L$(BOAT_DIR)/build $(LIBS)

clean:
	rm -f bo_naive_opti
