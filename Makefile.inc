# Build configuration.
#
# This file was originally set up for the author's Linux/cluster boxes:
#   LDFLAGS=-L/usr/lib -L/sw/lib -lgsl -lgslcblas
#   INCLUDES=-I/sw/include/
#   FLAGS=$(LDFLAGS) $(INCLUDES)
#
# Ported to build on macOS with GCC (the code is C++98-era: it uses
# <tr1/random>, dynamic exception specifications and dependent-base
# name lookup, so it is built with -std=gnu++98).
#
# GSL is located automatically: first in the bundled build under
# .deps/gsl (create it with 'make gsl'), otherwise in Homebrew's
# /opt/homebrew (Apple silicon) or /usr/local (Intel).
# See PORTING-NOTES.md for the full story of the macOS port.

DIRS=tools

# repository root, works both for 'make' here and 'make -C tools'
ROOT:=$(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))

GSL_PREFIX?=$(ROOT)/.deps/gsl
ifeq ($(wildcard $(GSL_PREFIX)/include/gsl),)
  GSL_INC?=-I/opt/homebrew/include
  GSL_LIB?=-L/opt/homebrew/lib
else
  GSL_INC=-I$(GSL_PREFIX)/include
  GSL_LIB=-L$(GSL_PREFIX)/lib
endif

# g++-14 (Homebrew GCC) still ships <tr1/random> and <tr1/unordered_map>;
# Apple clang / libc++ dropped TR1 entirely. Only override make's built-in
# default, so an explicit CC=... on the command line still wins.
ifeq ($(origin CC),default)
  CC:=$(shell command -v g++-14 2>/dev/null || command -v g++)
endif

CXXFLAGS_STD=-std=gnu++17
INCLUDES=$(GSL_INC)
FLAGS=$(CXXFLAGS_STD) $(INCLUDES)
LDFLAGS=$(GSL_LIB) -lgsl -lgslcblas -lm
LIBS=-L/usr/lib
SHELL=bash
ECHO=echo
RM=rm
MAKE=make

# Optimisation flags.  This used to be called DEBUGFLAGS and set to '-Wall -O2',
# which is neither debug nor particularly optimised.
#
# -O3 is the default.  -march=native is deliberately opt-in (NATIVE=1) because
# it makes the binary non-portable: it is fine for a benchmark on your own
# machine and wrong for anything you hand to somebody else.
WARNFLAGS=-Wall
OPTFLAGS?=-O3
NATIVE?=0
ifeq ($(NATIVE),1)
  OPTFLAGS+=-march=native
endif

# kept so that an existing 'make DEBUGFLAGS=...' still works
DEBUGFLAGS=$(WARNFLAGS) $(OPTFLAGS)
