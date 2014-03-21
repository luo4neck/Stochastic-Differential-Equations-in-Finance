// this file solves dX = aXdt + bXdW with both euler and milstein method.
// a and b are constants, I set a=b=2 here.

// X0 is red line with dt = 0.001 | milstein method
// X1 is green line with dt = 0.001 | euler method

// X2 is blue line with dt = 0.01 | milstein method
// X3 is purple line with dt = 0.01 | euler method

// X4 is purple line with dt = 0.1 | milstein method
// X5 is purple line with dt = 0.1 | euler method

#include<iostream>
#include<fstream>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_cdf.h>
#include<cmath>
#define SIZE 512
#define M 100

using namespace std;

class lines
{
    public:
    int size;
    double X[SIZE];
    double err[M];
    double mean;

    lines(int s): size(s)
    {
        mean = 0;
        for(int i=0; i<s; ++i)
        { X[i] = 0; }
        X[0] = 1;
        
        for(int i=0; i<M; ++i)
        { err[i] = 0; }
    }

    ~lines() {}

    double sigma()
    {
        double sigma=0;
        for(int i=0; i<M; ++i)
            sigma = sigma + (mean - err[i])*(mean - err[i]);
        
        return sigma/(double)(M-1); 
    }

    double up(double p)
    { return mean + gsl_cdf_tdist_Pinv(p, M-1)*sqrt( sigma()/(double)M); }
    
    double down(double p)
    { return mean + gsl_cdf_tdist_Qinv(p, M-1)*sqrt( sigma()/(double)M); }
};


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

int main()
{
gsl_rng *r = gsl_rng_alloc( gsl_rng_mt19937 );
gsl_rng_env_setup();
gsl_rng_set(r, time(NULL));
ofstream error("error.dat"); 
lines mil3(SIZE), elr3(SIZE);

double X0[SIZE], X1[SIZE], a=2, b=2, dt = 1.0/(double) SIZE, Xtrue[SIZE], W[SIZE], dW[SIZE];

for(int j=0; j<M; ++j)
{
Xtrue[0] = X0[0] = X1[0] = 1; 
mil3.err[j] = elr3.err[j] = 0;
Prepare(r, dt, a, b, dW, W, Xtrue);    
    for(int i=1; i<SIZE; ++i)
    {
        X0[i]= X0[i-1] + a*X0[i-1]*dt + b*X0[i-1]*dW[i-1] + 0.5*b*X0[i-1]*b*( dW[i-1]*dW[i-1] - dt);// milstein method..
        mil3.err[j] = mil3.err[j] + abs(Xtrue[i] - X0[i]);// error at this step..
        
        X1[i]= X1[i-1] + a*X1[i-1]*dt + b*X1[i-1]*dW[i-1]; // euler-maruyama method.. 
        elr3.err[j] = elr3.err[j] + abs(Xtrue[i] - X1[i]); // error at this step..
     /*   if((i+1)%10 == 0)
        {
            if((i+1)%100 == 0)
            {
            file<<(double)i/SIZE<<" "<<X1[i]<<" "<<X1[i]<<" "<<X1[i]<<" "<<0<<endl;
            }
            else 
            file<<(double)i/SIZE<<" "<<X1[i]<<" "<<X1[i]<<endl;
        }
        else*/
    }
mil3.err[j] = mil3.err[j]/(double)SIZE;
elr3.err[j] = elr3.err[j]/(double)SIZE;
mil3.mean = mil3.mean + mil3.err[j]; //accumulate the mean..
elr3.mean = elr3.mean + elr3.err[j];
error<<j<<" "<<mil3.err[j]<<" "<<elr3.err[j]<<endl;// print the error of two methods via gnuplot..
//cout<<sumXtrue<<endl<<sumX0<<endl<<sumX1<<endl;
}
gsl_rng_free(r);
error.close();
mil3.mean = mil3.mean/(double)M;//get the meam..
elr3.mean = elr3.mean/(double)M;


cout<<"The 90% confidence interval of milstein is:"<<endl;
cout<<mil3.up(0.95)<<" : "<<mil3.down(0.95)<<endl;

cout<<"The 90% confidence interval of euler is:"<<endl;
cout<<elr3.up(0.95)<<" : "<<elr3.down(0.95)<<endl;
/* ========== data store part ended, plotting part start ========== */

FILE *gp = popen("gnuplot -persist", "w");

if(gp == NULL)
{
cout<<"Cannot plot the data!"<<endl;
exit(0);
}

//fprintf(gp, "set logscale xy\n");
//fprintf(gp, "plot 'error.dat' u 1:2 w l, 'error.dat' u 1:3 w l, 'error.dat' u 1:4 w l, 'error.dat' u 1:5 w l\n");
//fprintf(gp, "plot 'error.dat' u 1:2 w l, 'error.dat' u 1:3 w l, 'error.dat' u 1:4 w l\n");
fprintf(gp, "plot 'error.dat' u 1:2 w l, 'error.dat' u 1:3 w l\n");
fprintf(gp, "pause -1\n");
fclose(gp);

//one side 90 95 97.5
//two side 80 90 95
return 0;
}
