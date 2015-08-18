
CXX=clang++

.PHONY: all clean

all: backprop_simple

backprop_simple: backprop_simple.cpp
	$(CXX) -O3 -std=c++14 -Wall -Werror -o $@ $<

clean:
	/bin/rm -rf backprop_simple

