
CC=g++

run: a.out
	./a.out

a.out:	*.cc *h ../ellipsoid.h
	$(CC)  main.cc -Wall  -I/sw/include/ -L/sw/lib/ -lgsl -lgslcblas


clean:
	rm -f a.out *.o
