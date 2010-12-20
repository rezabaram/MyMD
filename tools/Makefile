include ../Makefile.inc
CC=g++
INCLUDES=-I/sw/include/ -I$(HOME)/workstation/mysrc/
LIBS=-L/sw/lib/ -L$(HOME)/workstation/mysrc/
FLAGS=-Wall -lmylibs -lgsl -lgslcblas

targetdir=../bin/

all: coord_convert packing_density fabric 

packing_density:	packing_density.cc  ../include/ellipsoid.h
	$(ECHO) here
	$(CC) -O2 -o $(targetdir)/packing_density $(INCLUDES) $(LIBS) $(FLAGS) packing_density.cc 

coord_convert: coord_convert.cc ../include/packing.h
		g++ -o $(targetdir)/coord_convert $(INCLUDES) $(LIBS) $(FLAGS) coord_convert.cc 

#asphericity:asphericity.cc
#	$(CC)  $(INCLUDES) $(LIBS) $(FLAGS) asphericity.cc -o asphericity

fabric: fabric.cc ../include/packing.h ../include/ellips_contact.h ../include/interaction.h ../include/contact_network.h
	$(CC)  $(INCLUDES) $(LIBS) $(FLAGS) fabric.cc -o $(targetdir)/fabric

clean:
	rm -f fabric coord_convert packing_density *.o

network:fabric
	$(targetdir)/fabric packing.dat
plot: network
	$(targetdir)/coord2pr3d  -x0 -y0 -z0  -s1 network | render -size 1000x1000 -jpeg > _out.jpg && display _out.jpg