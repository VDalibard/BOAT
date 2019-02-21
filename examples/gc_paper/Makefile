CC=g++
CFLAGS=-Wall -fexceptions -std=c++14 -O3 -fopenmp
BOAT_DIR=./../../
LIBS= -lboat -lpthread -lblas -lnlopt -lm -lboost_serialization -lboost_system -lboost_coroutine -I $(EIGEN_DIR)

gc: gc.cpp
	$(CC) -o gc $(CFLAGS) -I$(BOAT_DIR)/include/ gc.cpp -L$(BOAT_DIR)/build $(LIBS)

clean:
	rm -f gc
