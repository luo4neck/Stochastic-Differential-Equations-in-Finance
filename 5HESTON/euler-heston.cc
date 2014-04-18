#include<iostream>
#include<fstream>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<cmath>
#define N 512 
#define MC 10000 // times of monte carlo tests..

using namespace std;

bool cor_check(double dw1[N], double dw2[N]) // check if the two processes are correlated..
{
    double sum12=0, sqsum1=0, sqsum2=0, sum1=0, sum2=0, w1[N], w2[N];
    sum1 = w1[0] = dw1[0], sum2 = w2[0] = dw2[0];
    sqsum1 = w1[0]*w1[0], sqsum2 = w2[0]*w2[0];
    sum12 = w1[0]*w2[0];

    for(int i=1; i<N; ++i)
    {
        w1[i] = w1[i-1] + dw1[i-1];
        sum1 = sum1 + w1[i];
        sqsum1 = sqsum1 + w1[i]*w1[i];

        w2[i] = w2[i-1] + dw2[i-1];
        sum2 = sum2 + w2[i];
        sqsum2 = sqsum2 + w2[i]*w2[i];
        
        sum12 = sum12 + w1[i]*w2[i];
    }
    double rho = (sum12/N - sum1*sum2/N/N) / sqrt( sqsum1/N - sum1*sum1/N/N ) / sqrt( sqsum2/N - sum2*sum2/N/N );
    
    if (rho > 0.4 ) 
    {
        return 0; // if rho is larger than 0.4, this path is usable, return 0 to exit the loop..
    }
    else return 1; // else return 1, generate dw1 and dw2 again and check the if correlated..
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
double strike_sum = 0;

for(int j=0; j<MC; ++j) // monte carlo part..
{
    bool check = 0;// check will decide to do correlation or not.. 0 is not do correlation.. 1 is do correlation..
    
    if (check == 0) generate_W(r, dW1, dW2, dt);
    
    while(check)
    {
        generate_W(r, dW1, dW2, dt);
        check = cor_check(dW1, dW2);
    }
    
    
    V[0] = 0.04, X[0] = 100;
    for(int i=1; i<N; ++i)
    {
        X[i] = X[i-1] + rate * X[i-1] * dt + sqrt(abs(V[i-1])) * X[i-1] * dW1[i-1];
        
        V[i] = V[i-1] + kappa * dt * (theta - V[i-1]) + eps * sqrt(abs(V[i-1])) * dW2[i-1];
    }
    
    if( X[N-1] > strike_price )
    {
        //strike_sum = strike_sum + (X[N-1] - strike_price) / pow( 1+rate , T);
        strike_sum = strike_sum + exp(-rate * T) * (X[N-1] - strike_price);
    }
}
gsl_rng_free(r);

cout<<"This is a Call Option at "<<strike_price<<"$"<<endl<<"Option price is "<<strike_sum/MC<<endl;

return 0;
}
