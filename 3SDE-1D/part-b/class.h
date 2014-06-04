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
