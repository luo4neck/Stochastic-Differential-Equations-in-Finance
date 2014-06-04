euler2D: euler2D.cc
	g++ -Wall euler2D.cc -O2 -o euler2D -lgsl -lgslcblas

test: euler2D 
	./euler2D

clean:
	rm -f *.dat core* euler2D milstein2D
