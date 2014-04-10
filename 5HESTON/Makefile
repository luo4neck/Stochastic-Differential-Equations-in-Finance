euler-heston: euler-heston.cc 
	g++ -Wall euler-heston.cc -O2 -o euler-heston -lgsl -lgslcblas

test: euler-heston
	./euler-heston

clean:
	rm -f *.o core* 
