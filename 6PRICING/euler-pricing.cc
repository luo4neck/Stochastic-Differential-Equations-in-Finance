#include<iostream>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<cmath>
#define N 512 
#define MC 10000 // times of monte carlo tests..

using namespace std;

void generate_W(gsl_rng *r, double w1[N], double w2[N], double dt)
{
    for(int i=0; i<N; ++i)
    {
       w1[i] = sqrt(dt) * gsl_ran_gaussian(r, 1.0); 
       w2[i] = sqrt(dt) * gsl_ran_gaussian(r, 1.0); 
    }
}

int main()
{
gsl_rng *r = gsl_rng_alloc( gsl_rng_mt19937);
gsl_rng_env_setup();
gsl_rng_set(r, time(NULL));

const double T = 1.0, rate = 0, dt = T/N;
double dW1[N], dW2[N], X[N], V[N]; 

const double strike_price = 1;
double strike_sum = 0; 

for(int j=0; j<MC; ++j) // monte carlo part..
{
    generate_W(r, dW1, dW2, dt); // generate two brownian motion paths..
    
    V[0] = 0.01, X[0] = 1;
    for(int i=1; i<N; ++i)// approach for euler SDE solution..
    {
        X[i] = X[i-1] + sqrt(max(V[i-1], 0.0)) * dW1[i-1]; // the above line shows that when V[i-1] is negative, set to 0, so right part at right of equal will be 0..
        //X[i] = X[i-1] + sqrt(abs(V[i-1])) * dW1[i-1]; // this line shows that take the absolute value of V[i-1]..
        
        V[i] = V[i-1] * ( 1 + dt - 0.5 * dW1[i-1] + 0.5 * sqrt(3) * dW2[i-1] ); 
    }
    
    if( X[N-1] > strike_price )
    {
        strike_sum = strike_sum + exp(-rate * T) * (X[N-1] - strike_price);
    }
}
gsl_rng_free(r);

cout<<"This is a Call Option at "<<strike_price<<"$"<<endl<<"Option price is "<<strike_sum/MC<<endl;

return 0;
}
