#Simple makefile for test purposes

CXX=g++
CXXFLAGS= -Wall -O2 -std=c++0x -g
LIBS= -lm
OUT=test_main
OBJS=vi_flux

$(OUT): $(OUT).cpp $(OBJS).o
	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@

$(OBJS): $(OBJS).cpp
	$(CXX) -c $(CXXFLAGS) $(LIBS)


clean:
	rm -f *.o *.dat *.png test_main

run:
	./main
