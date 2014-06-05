bpath: bpath.o 
	g++ -o bpath bpath.o -O2 -lgsl -lgslcblas 

bpath.o: bpath.cc 
	g++ -c bpath.cc -O2 -lgsl -lgslcblas

test: bpath
	./bpath

clean:
	rm -f *.o *.dat bpath 
