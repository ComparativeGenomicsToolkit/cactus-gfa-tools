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

all: mzgaf2paf pafcoverage rgfa-split paf2lastz rgfa2paf pafmask

mzgaf2paf: mzgaf2paf.o mzgaf2paf_main.o pafcoverage.o 
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -o mzgaf2paf mzgaf2paf_main.o mzgaf2paf.o pafcoverage.o $(LIB_FLAGS)

mzgaf2paf_main.o:$(LIB_DEPS) mzgaf2paf_main.cpp mzgaf2paf.hpp mzgaf.hpp gafkluge.hpp
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c mzgaf2paf_main.cpp $(INC_FLAGS)

mzgaf2paf.o:$(LIB_DEPS) mzgaf2paf.cpp mzgaf2paf.hpp mzgaf.hpp gafkluge.hpp pafcoverage.hpp
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c mzgaf2paf.cpp $(INC_FLAGS)

pafcoverage: pafcoverage.o pafcoverage_main.o
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -o pafcoverage pafcoverage_main.o pafcoverage.o $(LIB_FLAGS)

pafcoverage_main.o:$(LIB_DEPS) pafcoverage_main.cpp pafcoverage.hpp
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c pafcoverage_main.cpp $(INC_FLAGS)

pafcoverage.o:$(LIB_DEPS) pafcoverage.cpp pafcoverage.hpp
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c pafcoverage.cpp $(INC_FLAGS)

rgfa-split: rgfa-split.o rgfa-split_main.o pafcoverage.o
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -o rgfa-split rgfa-split_main.o rgfa-split.o pafcoverage.o $(LIB_FLAGS)

rgfa-split_main.o:$(LIB_DEPS) rgfa-split_main.cpp rgfa-split.hpp
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c rgfa-split_main.cpp $(INC_FLAGS)

rgfa-split.o:$(LIB_DEPS) rgfa-split.cpp rgfa-split.hpp gafkluge.hpp gfakluge.hpp pliib.hpp tinyfa.hpp pafcoverage.hpp
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c rgfa-split.cpp $(INC_FLAGS)

paf2lastz: paf2lastz.o paf2lastz_main.o
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -o paf2lastz paf2lastz_main.o paf2lastz.o

paf2lastz_main.o:$(LIB_DEPS) paf2lastz_main.cpp paf2lastz.hpp paf2lastz.hpp
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c paf2lastz_main.cpp $(INC_FLAGS)

paf2lastz.o:$(LIB_DEPS) paf2lastz.cpp paf2lastz.hpp
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c paf2lastz.cpp $(INC_FLAGS)

rgfa2paf_main.o: rgfa2paf_main.cpp gfakluge.hpp pafcoverage.hpp
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c rgfa2paf_main.cpp $(INC_FLAGS)

rgfa2paf: rgfa2paf_main.o pafcoverage.o
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -o rgfa2paf rgfa2paf_main.o pafcoverage.o

pafmask: pafmask_main.cpp rgfa-split.o pafcoverage.o
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c pafmask_main.cpp $(INC_FLAGS)
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -o pafmask pafmask_main.o pafcoverage.o rgfa-split.o

test : all paf2lastz_test pafmask_test
	cd test && prove -v test.t

paf2lastz_test: mapqTest scoreTest
	rm -f test/paf2lastz/out_mapq test/paf2lastz/out_score

mapqTest:
	./paf2lastz test/paf2lastz/evolver_rat.paf -q > test/paf2lastz/out_mapq
	diff test/paf2lastz/out_mapq test/paf2lastz/evolver_rat_mapq.cigar
	echo "OK"

scoreTest:
	./paf2lastz test/paf2lastz/evolver_rat.paf > test/paf2lastz/out_score
	diff test/paf2lastz/out_score test/paf2lastz/evolver_rat_score.cigar
	echo "OK"

pafmask_test:
	cd test && prove -v pafmask.t

clean:
	rm -rf mzgaf2paf main.o mzgaf2paf.o pafcoverage pafcoverage.o pafcoverage_main.o rgfa-split rgfa-split.o rgfa-split_main.o paf2lastz paf2lastz_main.o paf2lastz.o pafmask pafmask_main.o
