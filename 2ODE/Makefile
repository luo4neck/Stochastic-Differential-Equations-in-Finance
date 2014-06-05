euler: euler.o 
	g++ -Wall -o euler euler.o -O2 -lgsl -lgslcblas 

euler.o: euler.cc 
	g++ -Wall -c euler.cc -O2 -lgsl -lgslcblas

test: euler
	./euler

clean:
	rm -f *.o core* plot.dat euler 
