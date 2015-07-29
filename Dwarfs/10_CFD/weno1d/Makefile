CXX=g++
CC=gcc
LD=g++
CXXFLAGS=-g -O3 -fopenmp -Wall

OS = $(shell uname)
ifeq ($(OS),Darwin)
	CXX=clang-omp++
	CC=clang-omp
	LD=clang-omp++
	CXXFLAGS=-g -O3 -fopenmp -Wall
endif
