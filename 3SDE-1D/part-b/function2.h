/*
3.b.ii part, the given function is:
dX(t) = a*(b-X(t))*dt + c*sqrt(X(t))*dW(t) 
Problem of this question is the X(t) may become negative, 
so the analytic solution should be complex values. 
Analysing the error and mean of complex value and plot it
is too complex so what I am doing here is make the negative
values of Xtrue[i] to 0. 
*/
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_cdf.h>
#include<cmath>
#define SIZE 512
#define M 50000

using namespace std;

const double a = 2.0, b = 1.5, c = 1.0;

double euler(double X, double dW, double dt) // euler-maruyama method iteration function..
{
    X = max(X,0.0);
    return X + a * (b - X) * dt + c * sqrt(X) * dW; 
}

double milstein(double X, double dW, double dt) // milstein method iteration function..
{   
    X = max(X,0.0);
    return X + a * (b - X) * dt + c * sqrt(X) * dW + 0.25 * c * c * (dW * dW - dt); 
}

void Prepare(gsl_rng *r, double dt, double dW[SIZE], double W[SIZE], double Xtrue[SIZE])
{// this function gives dW and W values and return the true solution of the SDE
    double sum = Xtrue[0] = 1.0, X;
    W[0] = dW[0] = sqrt(dt) * gsl_ran_gaussian(r, 0);
    for(int i=1; i<SIZE; ++i)
    {
        dW[i] = sqrt(dt) * gsl_ran_gaussian(r, 1.0);
        W[i] = W[i-1] + dW[i-1];
        
        X = max(Xtrue[i-1], 0.0);
        Xtrue[i] = milstein(X, dW[i-1], dt); 
        // the actual solution.. 
        sum = sum + Xtrue[i];
    }
}
