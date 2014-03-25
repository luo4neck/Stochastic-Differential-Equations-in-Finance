// this file solves dX = aXdt + bXdW with both euler and milstein method.
// a and b are constants, I set a=b=2 here.

// mil3.X is red line with dt = 1/512 | milstein method
// elr3.X is green line with dt = 1/512 | euler method

// mil2.X is blue line with dt = 1/256 | milstein method
// elr2.X is purple line with dt = 1/256 | euler method

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
    double wmean;

    Lines(int s): size(s)
    {
        X = new double[s];
        mean = 0;
        wmean = 0;
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
Lines mil3(SIZE), elr3(SIZE), mil2(SIZE/2), elr2(SIZE/2), mil1(SIZE/4), elr1(SIZE/4), mil0(SIZE/8), elr0(SIZE/8);

double a=2, b=2, sum=0, dt = 1.0/(double) SIZE, Xtrue[SIZE], W[SIZE], dW[SIZE];
int I;

for(int j=0; j<M; ++j)
{
elr3.wmean = elr2.wmean = elr1.wmean = elr0.wmean = Xtrue[0] =1;

Prepare(r, dt, a, b, dW, W, Xtrue);    
    for(int i=1; i<SIZE; ++i)
    { // step size 1/512..
        mil3.X[i]= mil3.X[i-1] + a*mil3.X[i-1]*dt + b*mil3.X[i-1]*dW[i-1] + 0.5*b*mil3.X[i-1]*b*( dW[i-1]*dW[i-1] - dt);
        mil3.err[j] = mil3.err[j] + abs(Xtrue[i] - mil3.X[i]);
        elr3.X[i]= elr3.X[i-1] + a*elr3.X[i-1]*dt + b*elr3.X[i-1]*dW[i-1];
        elr3.err[j] = elr3.err[j] + abs(Xtrue[i] - elr3.X[i]); 
        elr3.wmean = elr3.wmean + elr3.X[i];
        
        if(i%2==0)
        {// step size 1/256..
        I = i/2;
        sum = dW[i-1] + dW[i-2];
        mil2.X[I]= mil2.X[I-1] + a*mil2.X[I-1]*2*dt + b*mil2.X[I-1]*sum + 0.5*b*mil2.X[I-1]*b*( sum*sum -2*dt);
        mil2.err[j] = mil2.err[j] + abs(Xtrue[I*2] - mil2.X[I]);
        elr2.X[I]= elr2.X[I-1] + a*elr2.X[I-1]*2*dt + b*elr2.X[I-1]*sum; 
        elr2.err[j] = elr2.err[j] + abs(Xtrue[I*2] - elr2.X[I]);
        elr2.wmean = elr2.wmean + elr2.X[I];
        
            if(i%4==0)
            {// step size 1/128..
            I = i/4;
            sum = dW[i-1] + dW[i-2] + dW[i-3] + dW[i-4];
            mil1.X[I]= mil1.X[I-1] + a*mil1.X[I-1]*4*dt + b*mil1.X[I-1]*sum + 0.5*b*mil1.X[I-1]*b*( sum*sum -4*dt);
            mil1.err[j] = mil1.err[j] + abs(Xtrue[I*4] - mil1.X[I]);
            elr1.X[I]= elr1.X[I-1] + a*elr1.X[I-1]*4*dt + b*elr1.X[I-1]*sum;
            elr1.err[j] = elr1.err[j] + abs(Xtrue[I*4] - elr1.X[I]);
            elr1.wmean = elr1.wmean + elr1.X[I];
            
                if(i%8==0)
                {// step size 1/64..
                I = i/8;
                sum = dW[i-1] + dW[i-2] + dW[i-3] + dW[i-4] + dW[i-5] + dW[i-6] + dW[i-7] + dW[i-8];
                mil0.X[I]= mil0.X[I-1] + a*mil0.X[I-1]*8*dt + b*mil0.X[I-1]*sum + 0.5*b*mil0.X[I-1]*b*( sum*sum -8*dt);
                mil0.err[j] = mil0.err[j] + abs(Xtrue[I*8] - mil0.X[I]);
                elr0.X[I]= elr0.X[I-1] + a*elr0.X[I-1]*8*dt + b*elr0.X[I-1]*sum;
                elr0.err[j] = elr0.err[j] + abs(Xtrue[I*8] - elr0.X[I]);
                elr0.wmean = elr0.wmean + elr0.X[I];
                }
            }
        }
    }
mil3.err[j] = mil3.err[j]/(double)mil3.size;
mil2.err[j] = mil2.err[j]/(double)mil2.size;
mil1.err[j] = mil1.err[j]/(double)mil1.size;
mil0.err[j] = mil0.err[j]/(double)mil0.size;
elr3.err[j] = elr3.err[j]/(double)elr3.size;
elr2.err[j] = elr2.err[j]/(double)elr2.size;
elr1.err[j] = elr1.err[j]/(double)elr1.size;
elr0.err[j] = elr0.err[j]/(double)elr0.size;

mil3.mean = mil3.mean + mil3.err[j]; //accumulate the mean..
mil2.mean = mil2.mean + mil2.err[j]; 
mil1.mean = mil1.mean + mil1.err[j];
mil0.mean = mil0.mean + mil0.err[j];
elr3.mean = elr3.mean + elr3.err[j];
elr2.mean = elr2.mean + elr2.err[j];
elr1.mean = elr1.mean + elr1.err[j];
elr0.mean = elr0.mean + elr0.err[j];

error<<j<<" "<<mil3.err[j]<<" "<<elr3.err[j]<<" "<<mil2.err[j]<<" "<<elr2.err[j]<<" "<<mil1.err[j]<<" "<<elr1.err[j]<<" "<<mil0.err[j]<<" "<<elr0.err[j]<<endl;// print the error of two methods via gnuplot..
}
gsl_rng_free(r);
error.close();

mil3.mean = mil3.mean/(double)M;//get the meam..
mil2.mean = mil2.mean/(double)M;
mil1.mean = mil1.mean/(double)M;
mil0.mean = mil0.mean/(double)M;

elr3.wmean = abs(elr3.wmean/(double)elr3.size/(double)M - exp(a));// weak convergence..
elr2.wmean = abs(elr2.wmean/(double)elr2.size/(double)M - exp(a));
elr1.wmean = abs(elr1.wmean/(double)elr1.size/(double)M - exp(a));
elr0.wmean = abs(elr0.wmean/(double)elr0.size/(double)M - exp(a));
cout<<elr3.wmean<<endl;
cout<<elr2.wmean<<endl;
cout<<elr1.wmean<<endl;
cout<<elr0.wmean<<endl;

elr3.mean = elr3.mean/(double)M;
elr2.mean = elr2.mean/(double)M;
elr1.mean = elr1.mean/(double)M;
elr0.mean = elr0.mean/(double)M;

ofstream data("data.dat");// send all data into 'data.dat' for error convergence..
data<<dt<<" "<<dt  <<" "<<mil3.mean<<" "<<mil3.up(0.975)<<" "<<mil3.down(0.975)<<" "<<mil3.up(0.95)<<" "<<mil3.down(0.95)<<" "<<elr3.mean<<" "<<elr3.up(0.975)<<" "<<elr3.down(0.975)<<" "<<elr3.up(0.95)<<" "<<elr3.down(0.95)<<endl;

data<<dt*2<<" "<<dt*2<<" "<<mil2.mean<<" "<<mil2.up(0.975)<<" "<<mil2.down(0.975)<<" "<<mil2.up(0.95)<<" "<<mil2.down(0.95)<<" "<<elr2.mean<<" "<<elr2.up(0.975)<<" "<<elr2.down(0.975)<<" "<<elr2.up(0.95)<<" "<<elr2.down(0.95)<<endl;

data<<dt*4<<" "<<dt*4<<" "<<mil1.mean<<" "<<mil1.up(0.975)<<" "<<mil1.down(0.975)<<" "<<mil1.up(0.95)<<" "<<mil1.down(0.95)<<" "<<elr1.mean<<" "<<elr1.up(0.975)<<" "<<elr1.down(0.975)<<" "<<elr1.up(0.95)<<" "<<elr1.down(0.95)<<endl;

data<<dt*8<<" "<<dt*8<<" "<<mil0.mean<<" "<<mil0.up(0.975)<<" "<<mil0.down(0.975)<<" "<<mil0.up(0.95)<<" "<<mil0.down(0.95)<<" "<<elr0.mean<<" "<<elr0.up(0.975)<<" "<<elr0.down(0.975)<<" "<<elr0.up(0.95)<<" "<<elr0.down(0.95)<<endl;
data.close();

/* ========== data store part ended, error plot part start ========== */

FILE *gp = popen("gnuplot -persist", "w");
if(gp == NULL)
{
cout<<"Cannot plot the data!"<<endl;
exit(0);
}
fprintf(gp, "set title '8 Lines show the error of 4 step sizes and 2 method'\n");
fprintf(gp, "plot 'error.dat' u 1:2 w l, 'error.dat' u 1:3 w l, 'error.dat' u 1:4 w l, 'error.dat' u 1:5 w l, 'error.dat' u 1:6 w l, 'error.dat' u 1:7 w l, 'error.dat' u 1:8 w l, 'error.dat' u 1:9 w l\n");
fprintf(gp, "pause -1\n");
fclose(gp);

/* ========== error plot ended, convergence plot start ========== */

FILE *fp = popen("gnuplot -persist", "w");
if(fp == NULL)
{
cout<<"Cannot plot the data!"<<endl;
exit(0);
}
fprintf(fp, "set logscale xy\n");
fprintf(fp, "set xrange [0.001:0.05]\n");

fprintf(fp, "f1(x)=a1*exp(x)+b1\n"); // fitting with exponential pattern..
fprintf(fp, "a1=7;b1=-7\n");
fprintf(fp, "f2(x)=a2*exp(x)+b2\n");
fprintf(fp, "a2=25;b2=-25\n");
/*
fprintf(fp, "f1(x)=a1*x+b1\n");// fitting with linear pattern..
fprintf(fp, "a1=7,b1=5\n");
fprintf(fp, "f2(x)=a2*x+b2\n");
fprintf(fp, "a2=5,b2=5\n");
*/
fprintf(fp, "fit [0.001:0.05] f1(x) 'data.dat' u 1:3 via a1,b1\n");
fprintf(fp, "fit [0.001:0.05] f2(x) 'data.dat' u 1:8 via a2,b2\n");
fprintf(fp, "plot 'data.dat' u 1:2 w l, 'data.dat' u 1:3:4:5 w yerrorlines, 'data.dat' u 1:3:6:7 w yerrorlines, 'data.dat' u 1:8:9:10 w yerrorlines, 'data.dat' u 1:8:11:12 w yerrorlines, f1(x) lw 1 lc rgb 'orange', f2(x) lw 1 lc rgb 'yellow'\n");

fprintf(fp, "pause -1\n");
fclose(fp); 
//one side 90 95 97.5
//two side 80 90 95
return 0;
}
