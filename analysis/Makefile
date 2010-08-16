
CC=g++

run: a.out
	./a.out

a.out:	*.cc *h ../ellipsoid.h
	$(CC)  main.cc -Wall  -I/sw/include/ -L/sw/lib/ -lgsl -lgslcblas

convert2raster: convert2raster.cc packing.h
		g++ -o convert2raster convert2raster.cc 

clean:
	rm -f a.out *.o
