
CC=g++
INCLUDES=-I/sw/include/ -I/Users/reza/workstation/mysrc/
LIBS=-L/sw/lib/ -L/Users/reza/workstation/mysrc/
FLAGS=-Wall -lmylibs -lgsl -lgslcblas

all: convert2raster packing_density 

packing_density:	*.cc *h ../ellipsoid.h
	$(CC) -o packing_density $(INCLUDES) $(LIBS) $(FLAGS) packing_density.cc 

convert2raster: convert2raster.cc packing.h
		g++ -o convert2raster convert2raster.cc 

asphericity:asphericity.cc
	$(CC)  $(INCLUDES) $(LIBS) $(FLAGS) asphericity.cc -o asphericity
clean:
	rm -f packing_density *.o
