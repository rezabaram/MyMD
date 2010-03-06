
run: a.out
	time ./run.sh

a.out:	*.cc *h
	g++ main.cc 

animate: test.avi
	mplayer test.avi

test.avi: 
	sh genFrames.sh out* > /dev/null 2>&1

clean:
	rm -rf a.out *o out* test.avi

