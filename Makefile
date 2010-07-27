
CC=g++

run: a.out
	time ./run.sh

a.out:	*.cc *h
	$(CC) -O2  main.cc -Wall  -I/sw/include/ -L/sw/lib/ -lgsl -lgslcblas

test:	
	plot.sh 1:2 log_energy trash/log_energy 
	#$(CC)  test.cc  -I/sw/include/ -L/sw/lib/ -lgsl -lgslcblas
diff:
	diff log_energy trash/log_energy
perf:
	plot.sh 1:2 log trash/log

animate: test.avi
	mplayer test.avi
	#feh *jpg

test.avi: 
	sh genFrames.sh out* > /dev/null 2>&1

clean:
	rm -f log_energy log out* *jpg test.avi 

mv:
	mv -f trash/* trash2/ &> /dev/null
	mv -f log_energy log out* *jpg test.avi trash &> /dev/null

aclean:
	rm -rf test.avi 

pov:
	coord2pov2  -s1.0 test.dat > test.pov && povray +h700 +w930 test.pov && feh test.png
	#coord2pov2  -s1.0 test.dat > test.pov && povray +h700 +w930 test.pov && feh test.png

zip:
	zip md.zip *.cc *h Makefile genFrames.sh run.sh config
