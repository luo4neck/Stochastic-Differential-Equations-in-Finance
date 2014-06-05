/*
the given function is:
dX(t) = aX(t)dt + bX(t)dW(t)
same as 3.a part
*/
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_cdf.h>
#include<cmath>
//#include "class.h"
#define SIZE 512
#define M 50000

using namespace std;

const double a = 2.0, b = 2.0;

double euler(double X, double dW, double dt) // euler-maruyama method iteration function..
{
    return (X + a * X * dt + b * X * dW);
}

double milstein(double X, double dW, double dt) // milstein method iteration function..
{
    return (X + a * X * dt + b * X * dW + 0.5 * b * X * b * (dW * dW - dt));
}

void Prepare(gsl_rng *r, double dt, double dW[SIZE], double W[SIZE], double Xtrue[SIZE])
//void Prepare(gsl_rng *r, double dt, double a, double b, double dW[SIZE], double W[SIZE], double Xtrue[SIZE])
{// this function gives dW and W values and return the true solution of the SDE
    double sum = Xtrue[0] = 1.0;
    W[0] = dW[0] = sqrt(dt) * gsl_ran_gaussian(r, 0);
    for(int i=1; i<SIZE; ++i)
    {
        dW[i] = sqrt(dt) * gsl_ran_gaussian(r, 1.0);
        W[i] = W[i-1] + dW[i-1];
        Xtrue[i] = Xtrue[0] * exp( (a-0.5*b*b) + b*W[i]); // the actual solution is a key part..
        sum = sum + Xtrue[i];
    }
}
