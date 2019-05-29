# IDIR is the directory for header files
IDIR =/usr/local/include
# g++ -L/usr/local/lib bessel_prova.o -lastrocpp -lgsl -lgslcblas -lm
# g++ -Wall -I/usr/local/include -c bessel_prova.cpp
# This is the compiler
CC=g++

# flags, i.e. option to pass to the compiler (tell where to search for stuff)
CFLAGS=-I$(IDIR)

# directory of libraries
LDIR =/usr/local/lib

# Option to put for compiling .o objects
CLIB=-L$(LDIR)

# directory of .o files
ODIR=obj

# libraries to add to the compiler
LIBS=-lm -lastrocpp -lgsl -lgslcblas -lKFT

# possible local header file
_DEPS =

# Put local header file in _DEPS and the one in include folder from IDIR, $ symbol is to call the variables
# you defined; % is for all the files under the IDIR directory;
# $(patsubst pattern,replacement,text) Finds whitespace-separated words in text that match pattern and
# replaces them with replacement; here every time you find something in _DEPS, you
# substitute with IDIR/something (so the output in this case will be DEPS = ../include/something.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

# The main cpp source that will be also the name of the program, usually change only this
MAIN_SOURCE = main

# put here also all the other sources
SOURCES = $(MAIN_SOURCE).cpp functions.cpp

_OBJ = $(patsubst %.cpp,%.o,$(SOURCES))
# Retrieving the object files specified in _OBJ in the subdirectory ODIR (here obj)
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

# targets are the object files in directory specified by ODIR, % is a special character (like * in bash);
# the prerequisites are the .c (.cpp,.cxx) files that have same name as the .o file (represented by %) and the
# header file in DEPS, that is compile looking at changes in the prerequisites
# $@ looks at the left of : whereas $< is the first prerequisite (here %.cpp, the source file)
$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

# $^ takes all prerequisites (the thing at the right of :)
$(MAIN_SOURCE): $(OBJ)
	$(CC) -o $@ $^ $(CLIB) $(LIBS)

.PHONY: clean
# clean the directory from all .o objects
clean:
	rm -f $(ODIR)/*.o
