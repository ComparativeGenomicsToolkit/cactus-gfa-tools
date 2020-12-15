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

CXXFLAGS := -O3 -Werror=return-type -std=c++14 -ggdb -g -MMD -MP $(PARALLEL_FLAGS) $(CXXFLAGS)

LIB_FLAGS = $(LIBS)
INC_FLAGS = -I$(CWD)

all: mzgaf2paf pafcoverage rgfa-split

mzgaf2paf: mzgaf2paf.o main.o
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -o mzgaf2paf main.o mzgaf2paf.o $(LIB_FLAGS)

main.o:$(LIB_DEPS) main.cpp mzgaf2paf.hpp mzgaf.hpp gafkluge.hpp
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c main.cpp $(INC_FLAGS)

mzgaf2paf.o:$(LIB_DEPS) mzgaf2paf.cpp mzgaf2paf.hpp mzgaf.hpp gafkluge.hpp
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c mzgaf2paf.cpp $(INC_FLAGS)

pafcoverage: pafcoverage.o pafcoverage_main.o
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -o pafcoverage pafcoverage_main.o pafcoverage.o $(LIB_FLAGS)

pafcoverage_main.o:$(LIB_DEPS) pafcoverage_main.cpp pafcoverage.hpp
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c pafcoverage_main.cpp $(INC_FLAGS)

pafcoverage.o:$(LIB_DEPS) pafcoverage.cpp pafcoverage.hpp
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c pafcoverage.cpp $(INC_FLAGS)

rgfa-split: rgfa-split.o rgfa-split_main.o
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -o rgfa-split rgfa-split_main.o rgfa-split.o $(LIB_FLAGS)

rgfa-split_main.o:$(LIB_DEPS) rgfa-split_main.cpp rgfa-split.hpp
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c rgfa-split_main.cpp $(INC_FLAGS)

rgfa-split.o:$(LIB_DEPS) rgfa-split.cpp rgfa-split.hpp gafkluge.hpp gfakluge.hpp pliib.hpp tinyfa.hpp
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c rgfa-split.cpp $(INC_FLAGS)

test : mzgaf2paf
	cd test && prove -v test.t

clean:
	rm -rf mzgaf2paf main.o mzgaf2paf.o pafcoverage pafcoverage.o pafcoverage_main.o rgfa-split rgfa-split.o rgfa-split_main.o
