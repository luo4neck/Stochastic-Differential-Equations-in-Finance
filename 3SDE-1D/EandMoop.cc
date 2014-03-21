// this file solves dX = aXdt + bXdW with both euler and milstein method.
// a and b are constants, I set a=b=2 here.

// mil3.X is red line with dt = 1/512 | milstein method
// elr3.X is green line with dt = 1/512 | euler method

// mil2.X is blue line with dt = 1/256 | milstein method
// elr2.X is purple line with dt = 1/256 | euler method

// mil1.X is purple line with dt = 1/128 | milstein method
// elr1.X is purple line with dt = 1/128 | euler method

#include<iostream>
#include<fstream>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_cdf.h>
#include<cmath>
#define SIZE 512
#define M 100

using namespace std;

class Lines
{
    private:
    double sigma;
    public:
    const int size;
    double *X;
    double err[M];
    double mean;

    Lines(int s): size(s)
    {
        X = new double[s];
        mean = 0;
        for(int i=0; i<s; ++i)
        { X[i] = 0; }
        X[0] = 1;
        
        for(int i=0; i<M; ++i)
        { err[i] = 0; }
    }

    ~Lines() {delete[] X;}

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
Lines mil3(SIZE), elr3(SIZE), mil2(SIZE/2), elr2(SIZE/2), mil1(SIZE/4), elr1(SIZE/4);

double a=2, b=2, sum=0, dt = 1.0/(double) SIZE, Xtrue[SIZE], W[SIZE], dW[SIZE];
int I;

for(int j=0; j<M; ++j)
{
Xtrue[0] =1;
Prepare(r, dt, a, b, dW, W, Xtrue);    
    for(int i=1; i<SIZE; ++i)
    {
        mil3.X[i]= mil3.X[i-1] + a*mil3.X[i-1]*dt + b*mil3.X[i-1]*dW[i-1] + 0.5*b*mil3.X[i-1]*b*( dW[i-1]*dW[i-1] - dt);
        // milstein method..
        mil3.err[j] = mil3.err[j] + abs(Xtrue[i] - mil3.X[i]);
        // error at this step..
        elr3.X[i]= elr3.X[i-1] + a*elr3.X[i-1]*dt + b*elr3.X[i-1]*dW[i-1];
        // euler-maruyama method.. 
        elr3.err[j] = elr3.err[j] + abs(Xtrue[i] - elr3.X[i]); 
        // error at this step..
        
        if(i%2==0)
        {
        I = i/2;
        sum = dW[i-1] + dW[i-2];
        mil2.X[I]= mil2.X[I-1] + a*mil2.X[I-1]*2*dt + b*mil2.X[I-1]*sum + 0.5*b*mil2.X[I-1]*b*( sum*sum -2*dt);
        mil2.err[j] = mil2.err[j] + abs(Xtrue[I*2] - mil2.X[I]);
        
        elr2.X[I]= elr2.X[I-1] + a*elr2.X[I-1]*2*dt + b*elr2.X[I-1]*sum; 
        elr2.err[j] = elr2.err[j] + abs(Xtrue[I*2] - elr2.X[I]);
        
            if(i%4==0)
            {
            I = i/4;
            sum = dW[i-1] + dW[i-2] + dW[i-3] + dW[i-4];
            mil1.X[I]= mil1.X[I-1] + a*mil1.X[I-1]*4*dt + b*mil1.X[I-1]*sum + 0.5*b*mil1.X[I-1]*b*( sum*sum -4*dt);
            mil1.err[j] = mil1.err[j] + abs(Xtrue[I*4] - mil1.X[I]);

            elr1.X[I]= elr1.X[I-1] + a*elr1.X[I-1]*4*dt + b*elr1.X[I-1]*sum;
            elr1.err[j] = elr1.err[j] + abs(Xtrue[I*4] - elr1.X[I]);
            //data<<I<<" "<<Xtrue[I*4]<<" "<<mil3.X[I*4]<<" "<<elr3.X[I*4]<<" "<<mil2.X[I*2]<<" "<<elr2.X[I*2]<<" "<<mil1.X[I]<<" "<<elr1.X[I]<<endl;// print the Xor of two methods via gnuplot..
            }
        }
    }
mil3.err[j] = mil3.err[j]/(double)mil3.size;
mil2.err[j] = mil2.err[j]/(double)mil2.size;
mil1.err[j] = mil1.err[j]/(double)mil1.size;
elr3.err[j] = elr3.err[j]/(double)elr3.size;
elr2.err[j] = elr2.err[j]/(double)elr2.size;
elr1.err[j] = elr1.err[j]/(double)elr1.size;

mil3.mean = mil3.mean + mil3.err[j]; //accumulate the mean..
mil2.mean = mil2.mean + mil2.err[j]; 
mil1.mean = mil1.mean + mil1.err[j];
elr3.mean = elr3.mean + elr3.err[j];
elr2.mean = elr2.mean + elr2.err[j];
elr1.mean = elr1.mean + elr1.err[j];

error<<j<<" "<<mil3.err[j]<<" "<<elr3.err[j]<<" "<<mil2.err[j]<<" "<<elr2.err[j]<<" "<<mil1.err[j]<<" "<<elr1.err[j]<<endl;// print the error of two methods via gnuplot..
}
gsl_rng_free(r);
error.close();

mil3.mean = mil3.mean/(double)M;//get the meam..
mil2.mean = mil2.mean/(double)M;
mil1.mean = mil1.mean/(double)M;
elr3.mean = elr3.mean/(double)M;
elr2.mean = elr2.mean/(double)M;
elr1.mean = elr1.mean/(double)M;

ofstream data("data.dat"); 
data<<1<<" "<<dt<<" "<<mil3.mean<<" "<<elr3.mean<<endl;
data<<2<<" "<<dt*2<<" "<<mil2.mean<<" "<<elr2.mean<<endl;
data<<3<<" "<<dt*4<<" "<<mil1.mean<<" "<<elr1.mean<<endl;
data.close();

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
fprintf(gp, "plot 'error.dat' u 1:2 w l, 'error.dat' u 1:3 w l, 'error.dat' u 1:4 w l, 'error.dat' u 1:5 w l, 'error.dat' u 1:6 w l, 'error.dat' u 1:7 w l\n");
//fprintf(gp, "plot 'error.dat' u 1:2 w l, 'error.dat' u 1:3 w l, 'error.dat' u 1:4 w l, 'error.dat' u 1:5 w l\n");
//fprintf(gp, "plot 'error.dat' u 1:2 w l, 'error.dat' u 1:3 w l, 'error.dat' u 1:4 w l\n");
//fprintf(gp, "plot 'error.dat' u 1:2 w l, 'error.dat' u 1:3 w l\n");
fprintf(gp, "pause -1\n");
fclose(gp);

FILE *fp = popen("gnuplot -persist", "w");
if(fp == NULL)
{
cout<<"Cannot plot the data!"<<endl;
exit(0);
}

fprintf(fp, "set logscale xy\n");
fprintf(fp, "plot 'data.dat' u 1:2 w l, 'data.dat' u 1:3 w l, 'data.dat' u 1:4 w l\n");
//fprintf(fp, "plot 'data.dat' u 1:2 w l, 'data.dat' u 1:3 w l\n");
//fprintf(fp, "plot 'data.dat' u 1:2 w l, 'data.dat' u 1:3 w l, 'data.dat' u 1:4 w l, 'data.dat' u 1:5 w l, 'data.dat' u 1:6 w l, 'data.dat' u 1:7 w l, 'data.dat' u 1:8 w l\n");
fprintf(fp, "pause -1\n");
fclose(fp); 
//one side 90 95 97.5
//two side 80 90 95
return 0;
}
