include Makefile.inc
DIRS= tools


run: ellipmd
	time bin/run.sh $(CONFIG)

ellipmd:	*.cc include/*.h 
	$(CC)  main.cc  $(FLAGS) -o ellipmd

condor_ellipmd:	*.cc include/*.h 
	$(CC)  main.cc /usr/lib64/gcc/x86_64-suse-linux/4.3/libstdc++.a  -Wall  $(LDFLAGS) $(DEBUGFLAGS) -o condor_ellipmd

tools:
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
	sh bin/genFrames.sh out0* > /dev/null 2>&1

movie:
	mencoder mf://*.jpg -mf fps=15:type=jpg -ovc lavc -ffourcc XVID -nosound -o out.avi

clean:
	rm -f log_energy log out0* *jpg test.avi 
	#$(ECHO) cleaning up in .
	#-$(RM) -f $(EXE) $(OBJS) $(OBJLIBS)
	#-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done


mv:
	mv -f trash/* trash2/ &> /dev/null
	mv -f log_energy log out* *jpg test.avi trash &> /dev/null

aclean:
	rm -rf test.avi 

pov:
	bin/coord_convert -p $(FILE) > out.dat && bin/coord2pov out.dat > out.pov && povray -A0.05 Antialias_Threshold=20  -w900 -h900 out.pov 

zip:
	zip md.zip *.cc *h Makefile genFrames.sh run.sh config

sync:
	rsync -ravz grace.cii.fc.ul.pt:workstation/MD/results/ results/

force_look :
	true
