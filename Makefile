
CC='g++'

run: a.out
	time ./run.sh

a.out:	*.cc *h
	$(CC)  main.cc -Wall  -I/sw64/include/ -L/sw64/lib/ -lgsl -lgslcblas

test:	test.cc 
	$(CC)  test.cc  -I/sw64/include/ -L/sw64/lib/ -lgsl -lgslcblas

animate: test.avi
	#feh *jpg
	mplayer test.avi

test.avi: 
	sh genFrames.sh out* > /dev/null 2>&1

clean:
	rm -rf a.out out* test.avi

aclean:
	rm -rf test.avi *jpg

pov:
	coord2pov2  -s1.0 test.dat > test.pov && povray +h700 +w930 test.pov && feh test.png
	#coord2pov2  -s1.0 test.dat > test.pov && povray +h700 +w930 test.pov && feh test.png

