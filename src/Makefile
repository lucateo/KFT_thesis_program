# IDIR is the directory for header files
IDIR =/usr/local/include
# This is the compiler
CC=g++

# flags, i.e. option to pass to the compiler (tell where to search for stuff)
CFLAGS=-I$(IDIR)

# directory of possible libraries
LDIR =/usr/local/lib

# Option to put for compiling .o objects
CLIB=-L$(LDIR)

# directory of .o files
ODIR=obj

# libraries to add to the compiler
LIBS=-lm -lastrocpp -lgsl -lgslcblas -lKFT

# possible local header file
DEPS =../include/KFT2.h

# Put local header file in _DEPS and the one in include folder from IDIR, $ symbol is to call the variables
# you defined; % is for all the files under the IDIR directory;
# $(patsubst pattern,replacement,text) Finds whitespace-separated words in text that match pattern and
# replaces them with replacement; here every time you find something in _DEPS, you
# substitute with IDIR/something (so the output in this case will be DEPS = ../include/something.h
#DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

# The main cpp source that will be also the name of the program, usually change only this
# Be aware of blank spaces!!!!!!!!!!!!!!!
# MAIN_SOURCE =mainTest
MAIN_SOURCE =mainCorrelation

# put here also all the other sources
SOURCES = $(MAIN_SOURCE).cpp KFT2.cpp functionTable.cpp osc_integrator.cpp bardeenPowerSpectrum.cpp smithEtalPowerSpectrum.cpp# functions.cpp write.cpp 

_OBJ = $(patsubst %.cpp,%.o,$(SOURCES))
# Retrieving the object files specified in _OBJ in the subdirectory ODIR (here obj)
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

# targets are the object files in directory specified by ODIR, % is a special character (like * in bash);
# the prerequisites are the .c (.cpp,.cxx) files that have same name as the .o file (represented by %) and the
# header file in DEPS, that is compile looking at changes in the prerequisites
# $@ looks at the left of : whereas $< is the first prerequisite (here %.cpp, the source file)
$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -g -o $@ $< $(CFLAGS)

# $^ takes all prerequisites (the thing at the right of :)
$(MAIN_SOURCE): $(OBJ)
	$(CC) -o $@ $^ $(CLIB) $(LIBS)

.PHONY: clean
# clean the directory from all .o objects
clean:
	rm -f $(ODIR)/*.o
