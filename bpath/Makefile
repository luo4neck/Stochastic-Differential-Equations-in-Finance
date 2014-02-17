bpath: bpath.o 
	g++ -o bpath bpath.o -lgsl -lgslcblas 

bpath.o: bpath.cc 
	g++ -c bpath.cc -lgsl -lgslcblas

test: bpath
	./bpath

clean:
	rm -f *.o 
