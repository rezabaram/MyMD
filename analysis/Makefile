
CC=g++
INCLUDES=-I/sw/include/ -I$(HOME)/workstation/mysrc/
LIBS=-L/sw/lib/ -L$(HOME)/workstation/mysrc/
FLAGS=-Wall -lmylibs -lgsl -lgslcblas

all: convert2raster packing_density fabric 

packing_density:	*.cc *h ../ellipsoid.h
	$(CC) -O2 -o packing_density $(INCLUDES) $(LIBS) $(FLAGS) packing_density.cc 

convert2raster: convert2raster.cc packing.h
		g++ -o convert2raster $(INCLUDES) $(LIBS) $(FLAGS) convert2raster.cc 

#asphericity:asphericity.cc
#	$(CC)  $(INCLUDES) $(LIBS) $(FLAGS) asphericity.cc -o asphericity

fabric: fabric.cc packing.h ../ellips_contact.h ../interaction.h contact_network.h
	$(CC)  $(INCLUDES) $(LIBS) $(FLAGS) fabric.cc -o fabric

clean:
	rm -f fabric convert2raster packing_density *.o

network:fabric
	./fabric packing.dat
plot: network
	coord2pr3d  -x0 -y0 -z0  -s1 network | render -size 1000x1000 -jpeg > _out.jpg && display _out.jpg
