all: euler-heston milstein-heston 

milstein-heston: milstein-heston.cc
	g++ -Wall milstein-heston.cc -O2 -o milstein-heston -lgsl -lgslcblas
	
euler-heston: euler-heston.cc
	g++ -Wall euler-heston.cc -O2 -o euler-heston -lgsl -lgslcblas

test: euler-heston milstein-heston
	./milstein-heston
	./euler-heston

clean:
	rm -f *.o core* euler-heston milstein-heston
