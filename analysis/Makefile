
CC=g++
INCLUDES=-I/sw/include/ -I/Users/reza/workstation/mysrc/
LIBS=-L/sw/lib/ -L/Users/reza/workstation/mysrc/
FLAGS=-Wall -lmylibs -lgsl -lgslcblas

run: a.out
	./a.out

a.out:	*.cc *h ../ellipsoid.h
	$(CC)  $(INCLUDES) $(LIBS) $(FLAGS) main.cc 

convert2raster: convert2raster.cc packing.h
		g++ -o convert2raster convert2raster.cc 

asphericity:asphericity.cc
	$(CC)  $(INCLUDES) $(LIBS) $(FLAGS) asphericity.cc -o asphericity
clean:
	rm -f a.out *.o
