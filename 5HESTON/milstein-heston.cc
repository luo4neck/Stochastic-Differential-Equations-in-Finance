#include<iostream>
#include<fstream>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<cmath>
#define N 512 
#define MC 10000 // times of monte carlo tests..

using namespace std;

void cor_check(double dw1[N], double dw2[N], double rho) // check if the two processes are correlated..
{
    for(int i=0; i<N; ++i)
        dw2[i] = dw1[i] * rho + dw2[i] * sqrt(1 - rho*rho);
}

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

const double T=1.0, rate = 0.05, kappa = 1.2, theta = 0.04, eps = 0.3, dt = T/N;
double dW1[N], dW2[N], X[N], V[N]; 

const double strike_price = 100;
const double rho=0.4; // this line set the correlation paramater rho = 0.4..
//const double rho=0.0;// this line set the correlation paramater rho = 0 ..
double strike_sum = 0; 

for(int j=0; j<MC; ++j) // monte carlo part..
{
    
    generate_W(r, dW1, dW2, dt); // generate two brownian paths..
    
    cor_check(dW1, dW2, rho); // correlation of the two brownian paths..
    
    V[0] = 0.04, X[0] = 100;
    for(int i=1; i<N; ++i) // milstein approach for SDE..
    {
        X[i] = X[i-1] + rate * X[i-1] * dt + sqrt(abs(V[i-1])) * X[i-1] * dW1[i-1] + 0.5 * X[i-1] * V[i-1] * (dW2[i-1]*dW2[i-1] - dt);
        
        V[i] = V[i-1] + kappa * dt * (theta - V[i-1]) + eps * sqrt(abs(V[i-1])) * dW2[i-1] + 0.25 * eps * eps * (dW2[i-1]*dW2[i-1] - dt); 
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
