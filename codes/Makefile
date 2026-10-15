ROOT_CONFIG ?= $(if $(ROOTSYS),$(ROOTSYS)/bin/root-config,root-config)
DELPHES_DIR ?= $(HOME)/wits/phd/software/Delphes-3.5.1

CXX ?= g++
CXXFLAGS ?= -O2 -Wall -Wextra
ROOTCFLAGS := $(shell $(ROOT_CONFIG) --cflags)
ROOTLIBS := $(shell $(ROOT_CONFIG) --libs)
CPPFLAGS += -I$(DELPHES_DIR)
CPPFLAGS += -DDELPHES_DIR='"$(DELPHES_DIR)"'
LDFLAGS += -L$(DELPHES_DIR) -Wl,-rpath,$(DELPHES_DIR)
# Delphes already contains ExRootAnalysis; do not link duplicate implementations.
LDLIBS += -lDelphes $(ROOTLIBS)

BASE ?= $(HOME)/wits/phd/ttaxion/bin
TARGET = $(BASE)/analyzedelphes
TEST_TARGET = $(BASE)/test_bjet_reference
SRC_HPP = utilities.hpp
SRC = utilities.cpp analyzedelphes.cpp

.PHONY: all dlph run test

all: $(TARGET)

dlph: all

$(TARGET): $(SRC) $(SRC_HPP) Makefile
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(ROOTCFLAGS) -std=c++17 $(SRC) $(LDFLAGS) $(LDLIBS) -o $@

run: all
	cd .. && $(abspath $(TARGET)) codes/params.json

$(TEST_TARGET): tests/test_bjet_reference.cpp utilities.cpp $(SRC_HPP) Makefile
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(ROOTCFLAGS) -std=c++17 utilities.cpp tests/test_bjet_reference.cpp $(LDFLAGS) $(LDLIBS) -o $@

test: $(TEST_TARGET)
	$(abspath $(TEST_TARGET))
