CXX=g++
CC=gcc
LD=g++
CXXFLAGS=-g -O3 -Wall -std=c++11
LDFLAGS=

OPENMP=-fopenmp

OS = $(shell uname)
ifeq ($(OS),Darwin)
	CXX=clang-omp++
	CC=clang-omp
	LD=clang-omp++
	CXXFLAGS=-g -O3 -Wall -std=c++11
	LDFLAGS=
endif

CXXFLAGS += $(OPENMP)
LDFLAGS += $(OPENMP)



all: solver1d

solver1d: solver1d.o
	$(LD) $(LDFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	/bin/rm -f solver1d solver1d.o
