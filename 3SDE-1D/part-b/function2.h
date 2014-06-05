/*
3.b.i part, the given function is:
dX(t) = a * X(t) * dt + b * dW(t)
This is a similar Langevin equation.
*/
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_cdf.h>
#include<cmath>
#define SIZE 512
#define M 50000

using namespace std;

const double a = 2.0, b = 0.5;

double euler(double X, double dW, double dt) // euler-maruyama method iteration function..
{
    return (X + a * X * dt + b * dW);
}

double milstein(double X, double dW, double dt) // milstein method iteration function..
{
    return (X + a * X * dt + b * dW); 
    // because the derivative of b is 0, the milstein part of this question is same with euler method..
}

void Prepare(gsl_rng *r, double dt, double dW[SIZE], double W[SIZE], double Xtrue[SIZE])
{// this function gives dW and W values and return the true solution of the SDE
    double sum = Xtrue[0] = 1.0, intgl = 0; 
    W[0] = dW[0] = sqrt(dt) * gsl_ran_gaussian(r, 0);
    for(int i=1; i<SIZE; ++i)
    {
        dW[i] = sqrt(dt) * gsl_ran_gaussian(r, 1.0);
        W[i] = W[i-1] + dW[i-1];
        intgl = intgl + exp(-1.0 * a * i * dt) * b * dW[i-1]; //calculate the integral part of actual solution.. 
        Xtrue[i] = exp(a * i * dt) * (Xtrue[0] + intgl); // the actual solution.. 
        sum = sum + Xtrue[i];
    }
}
