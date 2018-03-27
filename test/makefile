SRCS = $(wildcard ./*.cc)
EXES = $(SRCS:.cc=)
RUNS = $(SRCS:.cc=;)

CXX = g++
CXX_FLAG = -std=c++11 -I../include -I../src -O3 #-g #
LD_FLAGS = -lglog #-undefined dynamic_lookup 

.PHONY: all clean run
test_bayes: test_bayes.cc
	$(CXX) $(CXX_FLAG) $<  $(LD_FLAGS) -o $@ ../src/model/gp.cc ../src/optimizer/powell.cc ../src/optimizer/bayesopt.cc ../src/optimizer/ga.cc
	chmod +x $@
	./$@

shared: ../src/pybind11/py_gp.cc
	$(CXX) $(CXX_FLAG) $(LD_FLAGS) $^ ../src/model/gp.cc ../src/optimizer/powell.cc -undefined dynamic_lookup -shared -fPIC `python3 -m pybind11 --includes` -o gp`python3-config --extension-suffix`

test_gp: test_gp.cc
	$(CXX) $(CXX_FLAG) $<  $(LD_FLAGS) -o $@ ../src/model/gp.cc ../src/optimizer/powell.cc
	chmod +x $@
	./$@

test_powell: test_powell.cc
	$(CXX) $(CXX_FLAG) $<  $(LD_FLAGS) -o $@ ../src/optimizer/powell.cc
	chmod +x $@
	./$@

test_ga: test_ga.cc
	$(CXX) $(CXX_FLAG) $<  $(LD_FLAGS) -o $@ ../src/optimizer/ga.cc
	chmod +x $@
	./$@
	

re: clean all

all: $(EXES)

run: $(EXES)
	$(RUNS)


%: %.cc
	$(CXX) $(CXX_FLAG) $<  $(LD_FLAGS) -o $@ ../src/model/gp.cc
	chmod +x $@
	./$@

clean:
	-rm -rf $(EXES)



