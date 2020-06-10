# DnaParams is a header-only library (but associated Plumed ColVars are not).
# As a header library, there is no need to build anything to use it.
# This Makefile is for building and running unit tests.

.PHONY: test clean
 
LIB_DIR=lib
TEST_DIR=tests
CPPFLAGS=-std=c++17 -O3 -Wall

# Automatic dependency tracking, https://stackoverflow.com/a/30680399
CPPFLAGS += -MMD
-include $(wildcard *.d)

INCLUDES=\
	-I$(LIB_DIR)/\
	-I$(LIB_DIR)/eigen-3.3.7/\
	-I$(LIB_DIR)/autodiff-0.5.10/\
	-IDnaParams\
	-Iutils\
	-I$(TEST_DIR)
 
SRCS=$(wildcard $(TEST_DIR)/*.cpp)

OBJS=$(SRCS:.cpp=.o)
DEPS=$(SRCS:.cpp=.d)

1d29.pdb:
	wget https://files.rcsb.org/download/1d29.pdb


test: $(OBJS) 1d29.pdb
	g++ $(CPPFLAGS) $(INCLUDES) -o test $(OBJS)
	./test

%.o : %.cpp
	g++ $(CPPFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJS) $(DEPS)
