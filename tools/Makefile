include ../Makefile.inc
targetdir=../bin/

all: coord_convert packing_density fabric coord2xdr periodic2full compact correlation
	cp coord_convert packing_density fabric coord2xdr compact correlation $(targetdir)

packing_density:	packing_density.cc  ../include/packing.h ../include/ellipsoid.h
	$(CC) -O2 -o packing_density  packing_density.cc  $(FLAGS)

coord_convert: coord_convert.cc ../include/packing.h
		g++ -o coord_convert   coord_convert.cc  $(FLAGS)

coord2xdr: coord2xdr.cc
	g++ -o coord2xdr   coord2xdr.cc  $(FLAGS)

correlation: correlation.cc
	g++ -o correlation   correlation.cc  $(FLAGS)

compact: compact.cc
	g++ -o compact   compact.cc  $(FLAGS)

periodic2full: periodic2full.cc
	g++ -o periodic2full  periodic2full.cc   $(FLAGS)
#asphericity:asphericity.cc
#	$(CC)  $(INCLUDES) $(LIBS) $(FLAGS) asphericity.cc -o asphericity

fabric: fabric.cc ../include/packing.h ../include/ellips_contact.h ../include/interaction.h ../include/contact_network.h ../include/eigsys.h ../include/ellipsoid.h
	$(CC)    fabric.cc -o fabric  $(FLAGS)

clean:
	rm -f fabric coord_convert packing_density *.o

network:fabric
	$(targetdir)/fabric packing.dat
plot: network
	$(targetdir)/coord2pr3d  -x0 -y0 -z0  -s1 network | render -size 1000x1000 -jpeg > _out.jpg && display _out.jpg
