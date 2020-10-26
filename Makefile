CWD:=$(shell pwd)
CXX ?= g++

# Multithreading with OpenMP.
PARALLEL_FLAGS=-fopenmp -pthread
LIBS=-L$(CWD)

ifeq ($(shell uname -s),Darwin)
    # Our compiler might be clang that lacks -fopenmp support.
    # Sniff that
    ifeq ($(strip $(shell $(MY_CXX) -fopenmp /dev/null -o/dev/null 2>&1 | grep fopenmp | wc -l)), 1)
        # The compiler complained about fopenmp instead of its nonsense input file.
        # We need to use the hard way of getting OpenMP not bundled with the compiler.

        # The compiler only needs to do the preprocessing
        PARALLEL_FLAGS = -Xpreprocessor -fopenmp -pthread

        ifeq ($(shell if [ -d /opt/local/lib/libomp ];then echo 1;else echo 0;fi), 1)
            # Use /opt/local/lib/libomp if present, because Macports installs libomp there.
            # Brew is supposed to put it somewhere the compiler can find it by default.
            LIBS += -L/opt/local/lib/libomp
            # And we need to find the includes. Homebrew puts them in the normal place
            # but Macports hides them in "libomp"
            PARALLEL_FLAGS += -I/opt/local/include/libomp
        endif

        # We also need to link it
        LIBS += -lomp

    endif
endif

CXXFLAGS := -O0 -Werror=return-type -std=c++14 -ggdb -g -MMD -MP $(PARALLEL_FLAGS) $(CXXFLAGS)

LIB_FLAGS = $(LIBS)
INC_FLAGS = -I$(CWD)

all: mzgaf2paf

mzgaf2paf: mzgaf2paf.o main.o
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -o mzgaf2paf main.o mzgaf2paf.o $(LIB_FLAGS)

main.o:$(LIB_DEPS) main.cpp mzgaf2paf.hpp mzgaf.hpp gafkluge.hpp
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c main.cpp $(INC_FLAGS)

mzgaf2paf.o:$(LIB_DEPS) mzgaf2paf.cpp mzgaf2paf.hpp mzgaf.hpp gafkluge.hpp
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c mzgaf2paf.cpp $(INC_FLAGS)

test : mzgaf2paf
	cd test && prove -v hla.t

clean:
	rm -rf mzgaf2paf main.o mzgaf2paf.o
