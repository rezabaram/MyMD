include ../Makefile.inc
targetdir=../bin/

# NOTE (macOS port): $(FLAGS) carries -std=gnu++98 and the GSL include path,
# $(LDFLAGS) carries the GSL libraries. 'correlation' and 'coord2xdr'
# additionally require HDF5 and Sun-RPC/XDR (hdf5.h, rpc/xdr.h) which are not
# part of a stock macOS toolchain, so they are excluded from 'all' but kept
# as explicit targets.

all: coord_convert packing_density fabric compact periodic2full
	cp coord_convert packing_density fabric compact periodic2full $(targetdir)

packing_density:	packing_density.cc  ../include/packing.h ../include/ellipsoid.h
	$(CC) -O2 $(FLAGS) -o packing_density  packing_density.cc  $(LDFLAGS)

coord_convert: coord_convert.cc ../include/packing.h
	$(CC) $(FLAGS) -o coord_convert   coord_convert.cc  $(LDFLAGS)

coord2xdr: coord2xdr.cc
	$(CC) $(FLAGS) -o coord2xdr   coord2xdr.cc  $(LDFLAGS)

correlation: correlation.cc
	$(CC) $(FLAGS) -o correlation   correlation.cc  $(LDFLAGS)

compact: compact.cc
	$(CC) $(FLAGS) -o compact   compact.cc  $(LDFLAGS)

periodic2full: periodic2full.cc
	$(CC) $(FLAGS) -o periodic2full  periodic2full.cc   $(LDFLAGS)
#asphericity:asphericity.cc
#	$(CC)  $(INCLUDES) $(LIBS) $(FLAGS) asphericity.cc -o asphericity

fabric: fabric.cc ../include/packing.h ../include/ellips_contact.h ../include/interaction.h ../include/contact_network.h ../include/eigsys.h ../include/ellipsoid.h
	$(CC)    fabric.cc $(FLAGS) -o fabric  $(LDFLAGS)

clean:
	rm -f fabric coord_convert packing_density coord2xdr correlation compact periodic2full *.o

network:fabric
	$(targetdir)/fabric packing.dat
plot: network
	$(targetdir)/coord2pr3d  -x0 -y0 -z0  -s1 network | render -size 1000x1000 -jpeg > _out.jpg && display _out.jpg
