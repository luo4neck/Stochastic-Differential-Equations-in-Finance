em: em.o 
	g++ -Wall -o em em.o -O2 -lgsl -lgslcblas 

em.o: em.cc 
	g++ -Wall -c em.cc -O2 -lgsl -lgslcblas

test: em
	./em

clean:
	rm -f *.o *.dat  
