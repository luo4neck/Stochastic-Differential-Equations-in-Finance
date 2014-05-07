all: euler-pricing 

euler-pricing: euler-pricing.cc
	g++ -Wall euler-pricing.cc -O2 -o euler-pricing -lgsl -lgslcblas

test: euler-pricing 
	./euler-pricing

clean:
	rm -f *.o core* euler-pricing 
