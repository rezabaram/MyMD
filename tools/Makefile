include ../Makefile.inc
targetdir=../bin/

# Analysis utilities that read the snapshot format written by ellipmd.
#
# $(FLAGS) carries the C++ standard and the GSL include path; $(LDFLAGS)
# carries the GSL libraries.  The cluster/Fink paths this used to take over
# from the environment are gone -- see Makefile.inc.
#
# Two tools were removed here because they cannot be built outside the
# author's original setup: 'correlation' and 'coord2xdr' both need HDF5 and
# Sun RPC/XDR (hdf5.h, rpc/xdr.h), which are not part of a stock macOS
# toolchain.  'asphericity' and 'sphere_map' were already unbuildable (they
# include headers that no longer exist).

all: coord_convert packing_density fabric compact periodic2full
	cp coord_convert packing_density fabric compact periodic2full $(targetdir)

packing_density:	packing_density.cc  ../include/packing.h ../include/ellipsoid.h
	$(CC) -O2 $(FLAGS) -o packing_density  packing_density.cc  $(LDFLAGS)

coord_convert: coord_convert.cc ../include/packing.h
	$(CC) $(FLAGS) -o coord_convert   coord_convert.cc  $(LDFLAGS)

compact: compact.cc
	$(CC) $(FLAGS) -o compact   compact.cc  $(LDFLAGS)

periodic2full: periodic2full.cc
	$(CC) $(FLAGS) -o periodic2full  periodic2full.cc   $(LDFLAGS)

fabric: fabric.cc ../include/packing.h ../include/ellips_contact.h ../include/interaction.h ../include/contact_network.h ../include/eigsys.h ../include/ellipsoid.h
	$(CC)    fabric.cc $(FLAGS) -o fabric  $(LDFLAGS)

clean:
	rm -f fabric coord_convert packing_density compact periodic2full *.o
