emweak: emweak.o 
	g++ -Wall -o emweak emweak.o -O2 -lgsl -lgslcblas 

emweak.o: emweak.cc 
	g++ -Wall -c emweak.cc -O2 -lgsl -lgslcblas

test: emweak
	./emweak

clean:
	rm -f *.o core* 
