include ../Makefile.inc
CC=g++
INCLUDES=-I/sw/include/ -I$(HOME)/workstation/mysrc/
LIBS=-L/sw/lib/ -L$(HOME)/workstation/mysrc/
FLAGS=-Wall -lmylibs -lgsl -lgslcblas 
targetdir=../bin/

all: coord_convert packing_density fabric coord2xdr
	cp coord_convert packing_density fabric coord2xdr $(targetdir)

packing_density:	packing_density.cc  ../include/packing.h ../include/ellipsoid.h
	$(ECHO) here
	$(CC) -O2 -o packing_density $(INCLUDES) $(LIBS) $(FLAGS) packing_density.cc 

coord_convert: coord_convert.cc ../include/packing.h
		g++ -o coord_convert $(INCLUDES) $(LIBS) $(FLAGS) coord_convert.cc 

coord2xdr: coord2xdr.cc
	g++ -o coord2xdr $(INCLUDES) $(LIBS) $(FLAGS) coord2xdr.cc 

#asphericity:asphericity.cc
#	$(CC)  $(INCLUDES) $(LIBS) $(FLAGS) asphericity.cc -o asphericity

fabric: fabric.cc ../include/packing.h ../include/ellips_contact.h ../include/interaction.h ../include/contact_network.h ../include/eigsys.h ../include/ellipsoid.h
	$(CC)  $(INCLUDES) $(LIBS) $(FLAGS) fabric.cc -o fabric

clean:
	rm -f fabric coord_convert packing_density *.o

network:fabric
	$(targetdir)/fabric packing.dat
plot: network
	$(targetdir)/coord2pr3d  -x0 -y0 -z0  -s1 network | render -size 1000x1000 -jpeg > _out.jpg && display _out.jpg
