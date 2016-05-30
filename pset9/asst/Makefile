#
# Modify HALIDE_DIR to the appropriate path on your machine.
#
# Special instructions for Mac users
# ==================================
# You need to or must have installed libpng through Macports or Homebrew.
# Assuming that the installation succeeded, you should be able to run
#
# The brew command for installing libpng is
# brew install libpng
#
# libpng-config --I_opts
# libpng-config --L_opts
#
# Please add the output of the above commands to the following variables:
# PNG_INC
# PNG_LIB

MKDIR	:= mkdir -p
RM		:= rm -f
CP		:= cp -f
CXX		:= g++ -O3 -std=c++11 -I.

# TODO: changeme to your halide path
HALIDE_DIR ?= $(HOME)/halide
HALIDE_LIB := $(HALIDE_DIR)/bin/libHalide.a

INC  := $(wildcard  *.h)
SRC  := $(wildcard  tutorial*.cpp)

# TODO: (OSX users) add the configuration flags for libpng
PNG_INC := -I/usr/local/Cellar/libpng/1.6.17/include/libpng16
PNG_LIB := -L/usr/local/Cellar/libpng/1.6.17/lib

all: a9 
	mkdir -p Output

CXXFLAGS := $(PNG_INC) -I$(HALIDE_DIR)/include/ -I. -g -Wall
LDFLAGS  := $(PNG_LIB) -L$(HALIDE_DIR)/bin/     -lz -lpthread -ldl -lncurses -lpng

%: %.cpp $(HALIDE_LIB)
	$(CXX) $(CXXFLAGS) $< $(HALIDE_LIB) $(LDFLAGS) -o $@

a9: a9_main.cpp a9.cpp reference_implementations.cpp a9.h $(HALIDE_LIB)
	$(CXX) $(CXXFLAGS) a9_main.cpp a9.cpp reference_implementations.cpp timing.cpp $(HALIDE_LIB) $(LDFLAGS) -o $@

run: a9
	./a9

clean:
	$(RM) Output/* a9
	$(RM) -rf *.dSYM
	rm -rf Output
