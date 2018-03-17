SRCS = $(wildcard ./*.cc)
EXES = $(SRCS:.cc=)
RUNS = $(SRCS:.cc=;)

CXX = g++
CXX_FLAG = -std=c++11 -I../include -I../src -g
LD_FLAGS = -lglog

.PHONY: all clean run

all: run

re: clean all


run: $(EXES)
	$(RUNS)


%: %.cc
	$(CXX) $(CXX_FLAG) $<  $(LD_FLAGS) -o $@  ../src/optimizer/ga.cc
	chmod +x $@

clean:
	-rm -rf $(EXES)



