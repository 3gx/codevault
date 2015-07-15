CXX = g++
CC  = gcc
LD  = g++
F90 = gfortran

ifeq ($(OMP), 1)
	OMPFLAGS  = -fopenmp
endif
ifeq ($(OMPGLIBC), 1)
	OMPFLAGS += -D_GLIBCXX_PARALLEL
endif

OFLAGS = -O4 -g -Wall
ifeq ($(GDB), 1)
	OFLAGS = -O0 -g -Wall
endif
	
OFLAGS  += -funroll-all-loops
CXXFLAGS = $(OFLAGS) -Wstrict-aliasing=2 $(OMPFLAGS)

ifeq ($(SSE), 1)
	CXXFLAGS += -D__mySSE__ -msse4.1
endif
ifeq ($(AVX), 1)
	CXXFLAGS += -D__myAVX__ -mavx -mtune=corei7-avx -march=corei7-avx
endif
ifeq ($(SSE1), 1)
	CXXFLAGS += -D__mySSE1__ -msse4.1
endif
ifeq ($(RELEASE), 1)
	CXXFLAGS += -DNDEBUG
endif



LDFLAGS = 


SRCPATH = ./
SRC = octree.cpp
OBJ = $(SRC:%.cpp=%.o) 

PROG  = octree

RM = /bin/rm

all:	  $(PROG) $(PROG1) $(PROG2)


$(PROG): $(OBJ) 
	$(LD) $(LDFLAGS) $^ -o $@ $(OMPFLAGS)

%.o: $(SRCPATH)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean:
	/bin/rm -rf *.o $(PROG) 

clean_all:
	/bin/rm -rf *.o $(PROG) plummer.dat *~

$(OBJ): mytimer.h sse.h avx.h vector3.h boundary4.h peano.h primitives.h \
	octree.h Body.h Cell.h Boundaries.h Group.h \
	range_search_ptcl.h range_search_group.h \
	Multipole.h gforce_group.h







