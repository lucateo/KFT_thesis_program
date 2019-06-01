# KFT_thesis_program
C++ program for the analysis of the non linear evolution of power spectra using the Kinetic Field Theory framework

This is the program I am attempting to develop for my thesis on KFT;

It uses the libraries libastro and libKFT, so to use it you should install those libraries first;
I put also the Makefile I used to compile, it assumes you are in the src directory, that your
libraries are in /usr/local/include and that you have the obj/ directory (to store .o files necessary for
compilation). It also uses g++. You can easily change all this by checking the 
Makefile, then you compile by simply typing 'make' ( 'make clean' to remove all .o files).
The code is documented with Doxygen comment blocks, so you can extract a Doxygen documentation from it.
