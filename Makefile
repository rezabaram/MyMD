include Makefile.inc
DIRS= tools


run: a.out
	time bin/run.sh

a.out:	*.cc include/*.h 
	$(CC)  main.cc -Wall  $(LDFLAGS)

subsystem:
	$(MAKE) -C tools

test:	
	plot.sh 1:2 log_energy trash/log_energy 
	#$(CC)  test.cc  -I/sw/include/ -L/sw/lib/ -lgsl -lgslcblas
diff:
	diff log_energy trash/log_energy
perf:
	plot.sh 1:2 log trash/log

animate: test.avi
	#mplayer -loop 0 test.avi
	feh *jpg

test.avi: 
	sh bin/genFrames.sh out* > /dev/null 2>&1

clean:
	rm -f log_energy log out* *jpg test.avi 
	#$(ECHO) cleaning up in .
	#-$(RM) -f $(EXE) $(OBJS) $(OBJLIBS)
	#-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done


mv:
	mv -f trash/* trash2/ &> /dev/null
	mv -f log_energy log out* *jpg test.avi trash &> /dev/null

aclean:
	rm -rf test.avi 

pov:
	analysis/coord_convert -p $(FILE) > out.dat && bin/coord2pov -c out.dat > out.pov && povray -A0.05 Antialias_Threshold=20  -w500 -h500 out.pov && feh out.png

zip:
	zip md.zip *.cc *h Makefile genFrames.sh run.sh config

sync:
	rsync -ravz grace.cii.fc.ul.pt:workstation/MD/results/ results

force_look :
	true

