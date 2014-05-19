#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_cdf.h>
#include<cmath>
#define SIZE 512
#define M 50000

using namespace std;

class Line_strong
{
    private:
    double sigma;
    public:
    const int size;
    double *X;
    double err[M];
    double mean;

    Line_strong(int s): size(s)
    {
        X = new double[s];
        mean = 0;
        for(int i=0; i<s; ++i)
        { X[i] = 0; }
        X[0] = 1;
        
        for(int i=0; i<M; ++i)
        { err[i] = 0; }
    }

    ~Line_strong() 
    {
        delete[] X;
    }

    double getsigma()
    {
        double S=0;
        for(int i=0; i<M; ++i)
            S = S + (mean - err[i])*(mean - err[i]);
        return S/(double)(M-1); 
    }

    double up(double p)
    {
        sigma = getsigma(); 
        return mean + gsl_cdf_tdist_Pinv(p, M-1)*sqrt( sigma/(double)M); 
    }
    double down(double p)
    {
        sigma = getsigma(); 
        return mean + gsl_cdf_tdist_Qinv(p, M-1)*sqrt( sigma/(double)M); 
    }
};

class Line_weak
{
    public:
    const int size;
    double *X;
    double mean;

    Line_weak(int s): size(s)
    {
        X = new double[s];
        mean = 0;
        for(int i=0; i<s; ++i)
        { X[i] = 0; }
        X[0] = 1;
    }

    ~Line_weak() 
    {
        delete[] X;
    }
};

double sign(double sum)
{
    if (sum < 0) return -1.0;
    else if (sum > 0) return 1.0;
    else return 0.0;
}

void Prepare(gsl_rng *r, double dt, double a, double b, double dW[SIZE], double W[SIZE], double Xtrue[SIZE])
{// this function gives dW and W values and return the true solution of the SDE
    double sum=Xtrue[0];
    W[0] = dW[0] = sqrt(dt) * gsl_ran_gaussian(r, 0);
    for(int i=1; i<SIZE; ++i)
    {
        dW[i] = sqrt(dt) * gsl_ran_gaussian(r, 1.0);
        W[i] = W[i-1] + dW[i-1];
        Xtrue[i] = 1.0 * exp( (a-0.5*b*b) + b*W[i]);
        sum = sum + Xtrue[i];
    }
}
