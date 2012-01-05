DIRS=tools 
CC=g++
LDFLAGS=-L/usr/lib -L/sw/lib -lgsl -lgslcblas
INCLUDES=-I/sw/include/
FLAGS=$(LDFLAGS) $(INCLUDES)
LIBS=-L/usr/lib
SHELL=bash
ECHO=echo
RM=rm
MAKE=make
DEBUGFLAGS=-Wall -O2

