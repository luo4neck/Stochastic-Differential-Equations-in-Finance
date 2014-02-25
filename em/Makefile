emstrong: emstrong.o 
	g++ -Wall -o emstrong emstrong.o -O2 -lgsl -lgslcblas 

emstrong.o: emstrong.cc 
	g++ -Wall -c emstrong.cc -O2 -lgsl -lgslcblas

test: emstrong
	./emstrong

clean:
	rm -f *.o core* 
